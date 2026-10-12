################################################################################
# 02-prepare-lesion-tables.R
#
# Subset the SNV consensus + hotspot MAF to the genes named in
# input/subtype_lesion_rules.tsv and input/annotation_genes.tsv and write a small parquet file, so that
# 03-validate-subtypes.R does not need to re-read the multi-GB MAF.
#
# The MAF is a single gzip stream, so it is scanned once here (single
# threaded) rather than on every validator run.
#
# Usage: Rscript --vanilla 02-prepare-lesion-tables.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS")
source(file.path(module_dir, "utils", "duckdb-connect.R"))
scratch_dir <- file.path(root_dir, "scratch", "molecular-subtyping-CNS")
dir.create(scratch_dir, showWarnings = FALSE, recursive = TRUE)

maf_file <- file.path(root_dir, "data", "snv-consensus-plus-hotspots.maf.tsv.gz")
out_file <- file.path(scratch_dir, "snv_lesion_genes.parquet")

rules <- read_tsv(file.path(module_dir, "input", "subtype_lesion_rules.tsv"),
                  show_col_types = FALSE)

annotation_genes <- read_tsv(file.path(module_dir, "input", "annotation_genes.tsv"),
                             show_col_types = FALSE)

lgg_components <- read_tsv(file.path(module_dir, "input", "lgg_component_rules.tsv"),
                           show_col_types = FALSE)

genes <- rules %>%
  filter(lesion_type == "snv") %>%
  pull(gene) %>%
  c(annotation_genes$gene,
    lgg_components$gene[lgg_components$lesion_type == "snv"]) %>%
  str_split("\\|") %>%
  unlist() %>%
  unique()

gene_sql <- paste0("'", genes, "'", collapse = ", ")

con <- connect_duckdb()
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

# all_varchar and a short sample avoid type sniffing over the whole file
query <- sprintf(
  "COPY (
     SELECT Hugo_Symbol, HGVSp_Short, Variant_Classification, DOMAINS,
            Tumor_Sample_Barcode AS Kids_First_Biospecimen_ID
     FROM read_csv('%s', delim = '\t', header = true, all_varchar = true,
                   sample_size = 100, ignore_errors = true)
     WHERE Hugo_Symbol IN (%s)
   ) TO '%s' (FORMAT PARQUET)",
  maf_file, gene_sql, out_file
)

start <- Sys.time()
dbExecute(con, query)
message("Wrote ", out_file, " in ",
        round(difftime(Sys.time(), start, units = "mins"), 1), " min")

# 1p/19q arm-level values from the consensus GISTIC results (WGS only), used
# for the 1p/19q codeletion rule. Only the two rows needed are kept.
arm_zip <- file.path(root_dir, "data", "cnv-consensus-gistic.zip")
arm_tmp <- tempfile("gistic")
dir.create(arm_tmp)
unzip(arm_zip, files = "cnv-consensus-gistic/broad_values_by_arm.txt",
      exdir = arm_tmp, junkpaths = TRUE)

arm_df <- data.table::fread(file.path(arm_tmp, "broad_values_by_arm.txt"),
                            data.table = FALSE) %>%
  filter(`Chromosome Arm` %in% c("1p", "19q")) %>%
  column_to_rownames("Chromosome Arm") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Kids_First_Biospecimen_ID")

write_tsv(arm_df, file.path(scratch_dir, "arm_1p19q.tsv"))
message("Wrote 1p/19q arm values for ", nrow(arm_df), " biospecimens")
unlink(arm_tmp, recursive = TRUE)

# Expression of the genes in input/expression_annotations.tsv (e.g. EZHIP
# overexpression), subset from the collapsed TPM matrix
expr_genes <- read_tsv(file.path(module_dir, "input", "expression_annotations.tsv"),
                       show_col_types = FALSE)
tpm <- readRDS(file.path(root_dir, "data", "gene-expression-rsem-tpm-collapsed.rds"))
stopifnot(all(expr_genes$gene %in% rownames(tpm)))

# Overexpression is a z-score of log2(TPM + 1) within RNA library type, against
# PBTA tumor RNA-seq biospecimens. Library type matters: exome-capture
# libraries have much higher EZHIP TPM than poly-A or stranded ones.
rna_meta <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                     guess_max = 100000, show_col_types = FALSE) %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  transmute(Kids_First_Biospecimen_ID, cohort, sample_type,
            library_group = case_when(
              str_detect(RNA_library, "exome") ~ "exome capture",
              str_detect(RNA_library, "poly-A") ~ "poly-A",
              TRUE ~ "stranded"
            ))

expr_long <- tpm[expr_genes$gene, , drop = FALSE] %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "Kids_First_Biospecimen_ID", values_to = "tpm") %>%
  inner_join(rna_meta, by = "Kids_First_Biospecimen_ID") %>%
  mutate(log2_tpm = log2(tpm + 1))

reference <- expr_long %>%
  filter(cohort == "PBTA", sample_type == "Tumor") %>%
  group_by(gene, library_group) %>%
  summarise(ref_mean = mean(log2_tpm), ref_sd = sd(log2_tpm), n_ref = n(),
            .groups = "drop")

expr_long <- expr_long %>%
  left_join(reference, by = c("gene", "library_group")) %>%
  mutate(z = (log2_tpm - ref_mean) / ref_sd) %>%
  select(gene, Kids_First_Biospecimen_ID, library_group, tpm, z)

write_tsv(expr_long, file.path(scratch_dir, "expression_annotation_genes.tsv"))
message("Wrote expression for ", nrow(expr_genes), " gene(s) in ",
        n_distinct(expr_long$Kids_First_Biospecimen_ID), " biospecimens")

# Gene-level copy-number status for the genes of cnv rules (e.g. MYC and GLI2
# amplification in medulloblastoma), subset from the consensus gene-level file
cnv_genes <- rules %>%
  filter(lesion_type == "cnv") %>%
  pull(gene) %>%
  str_split("\\|") %>%
  unlist() %>%
  unique()

if (length(cnv_genes) > 0) {
  cnv_file <- file.path(root_dir, "data",
    "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_autosomes.tsv.gz")
  cnv_out <- file.path(scratch_dir, "cnv_lesion_genes.parquet")
  cnv_query <- sprintf(
    "COPY (
       SELECT biospecimen_id AS Kids_First_Biospecimen_ID, gene_symbol, status
       FROM read_csv('%s', delim = '\t', header = true, all_varchar = true,
                     sample_size = 100, ignore_errors = true)
       WHERE gene_symbol IN (%s) AND status NOT IN ('neutral', 'unchanged')
     ) TO '%s' (FORMAT PARQUET)",
    cnv_file, paste0("'", cnv_genes, "'", collapse = ", "), cnv_out
  )
  con2 <- connect_duckdb()
  dbExecute(con2, cnv_query)
  dbDisconnect(con2, shutdown = TRUE)
  message("Wrote gene-level copy number for ", length(cnv_genes), " gene(s)")
}

# LGG module lesions that need more than a gene list ------------------------------
# CDKN2A / CDKN2B copy number (legacy LGG rule: copy number <= 1)
cdkn_file <- file.path(root_dir, "data",
  "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz")
cdkn_out <- file.path(scratch_dir, "cdkn2_copy_number.parquet")
con3 <- connect_duckdb()
dbExecute(con3, sprintf(
  "COPY (
     SELECT biospecimen_id AS Kids_First_Biospecimen_ID, gene_symbol,
            TRY_CAST(copy_number AS DOUBLE) AS copy_number, status
     FROM read_csv('%s', delim = '\t', header = true, all_varchar = true,
                   sample_size = 100, ignore_errors = true)
     WHERE gene_symbol IN ('CDKN2A', 'CDKN2B')
   ) TO '%s' (FORMAT PARQUET)", cdkn_file, cdkn_out))

# FGFR1 kinase-domain duplication, as in the LGG module: a tandem duplication
# (Manta) that spans the FGFR1 protein tyrosine kinase domain, or a consensus
# copy-number segment with copy number > 2 within the FGFR1 gene
pfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS",
                            package = "annoFuseData")) %>%
  filter(hgnc_symbol == "FGFR1", DESC == "Protein tyrosine kinase",
         !is.na(domain_start), !is.na(domain_end)) %>%
  distinct(chromosome_name, domain_start, domain_end, gene_start, gene_end)
stopifnot(nrow(pfam) == 1)

manta_out <- file.path(scratch_dir, "fgfr1_tandem_dup.parquet")
dbExecute(con3, sprintf(
  "COPY (
     SELECT \"Kids.First.Biospecimen.ID.Tumor\" AS Kids_First_Biospecimen_ID,
            SV_start, SV_end, ALT
     FROM read_csv('%s', delim = '\t', header = true, all_varchar = true,
                   sample_size = 100, ignore_errors = true)
     WHERE Gene_name = 'FGFR1' AND ALT = '<DUP:TANDEM>'
   ) TO '%s' (FORMAT PARQUET)",
  file.path(root_dir, "data", "sv-manta.tsv.gz"), manta_out))

manta_fgfr <- dbGetQuery(con3, sprintf("SELECT * FROM read_parquet('%s')", manta_out)) %>%
  as_tibble() %>%
  mutate(across(c(SV_start, SV_end), as.numeric)) %>%
  filter(SV_start <= pfam$domain_start, SV_end >= pfam$domain_end) %>%
  distinct(Kids_First_Biospecimen_ID) %>%
  mutate(fgfr1_dup_source = "tandem duplication (SV)")

seg_fgfr <- read_tsv(file.path(root_dir, "data", "cnv-consensus.seg.gz"),
                     show_col_types = FALSE) %>%
  filter(chrom == paste0("chr", pfam$chromosome_name),
         loc.start >= pfam$gene_start, loc.end <= pfam$gene_end,
         copy.num > 2) %>%
  distinct(Kids_First_Biospecimen_ID = ID) %>%
  mutate(fgfr1_dup_source = "copy-number segment within FGFR1")

bind_rows(manta_fgfr, seg_fgfr) %>%
  group_by(Kids_First_Biospecimen_ID) %>%
  summarise(fgfr1_dup_source = paste(sort(unique(fgfr1_dup_source)), collapse = "; "),
            .groups = "drop") %>%
  write_tsv(file.path(scratch_dir, "fgfr1_dup.tsv"))
dbDisconnect(con3, shutdown = TRUE)
message("Wrote CDKN2A/B copy number and FGFR1 duplication calls")

# C19MC region copy number (embryonal tumors, ETMR C19MC-altered) --------------
# chr19 segments with a positive segment mean from the consensus copy number;
# 03 keeps the ones that overlap the C19MC miRNA cluster (hg38 chr19:53.66-53.77 Mb)
read_tsv(file.path(root_dir, "data", "cnv-consensus.seg.gz"),
         show_col_types = FALSE) %>%
  filter(chrom == "chr19", seg.mean > 0) %>%
  transmute(Kids_First_Biospecimen_ID = ID, loc.start, loc.end, seg.mean, copy.num) %>%
  write_tsv(file.path(scratch_dir, "chr19_gain_segments.tsv"))
message("Wrote chr19 gain segments")
