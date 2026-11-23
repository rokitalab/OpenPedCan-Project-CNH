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
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
source(file.path(module_dir, "utils", "duckdb-connect.R"))
source(file.path(module_dir, "utils", "release-tables.R"))
source(file.path(module_dir, "utils", "fusion-calls.R"))
suppressPackageStartupMessages(library(annoFuseData))
scratch_dir <- file.path(root_dir, "scratch", "molecular-subtyping")
dir.create(scratch_dir, showWarnings = FALSE, recursive = TRUE)

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

# the parquet copy of the MAF is sorted by gene, so only the matching row
# groups are read; the TSV fallback is read as text to avoid type sniffing
query <- sprintf(
  "COPY (
     SELECT Hugo_Symbol,
            -- TERT promoter hotspots have no protein change: named by the
            -- classic cDNA labels (chr5:1295113 G>A is C228T, 1295135 G>A C250T)
            CASE WHEN Hugo_Symbol = 'TERT' AND Chromosome = 'chr5' AND Tumor_Seq_Allele2 = 'A'
                      AND Start_Position = 1295113 THEN 'promoter C228T'
                 WHEN Hugo_Symbol = 'TERT' AND Chromosome = 'chr5' AND Tumor_Seq_Allele2 = 'A'
                      AND Start_Position = 1295135 THEN 'promoter C250T'
                 ELSE HGVSp_Short END AS HGVSp_Short,
            Variant_Classification, DOMAINS,
            Tumor_Sample_Barcode AS Kids_First_Biospecimen_ID
     FROM %s
     WHERE Hugo_Symbol IN (%s)
   ) TO '%s' (FORMAT PARQUET)",
  release_source("snv-consensus-plus-hotspots",
                 "snv-consensus-plus-hotspots.maf.tsv.gz", root_dir),
  gene_sql, out_file
)

start <- Sys.time()
dbExecute(con, query)
message("Wrote ", out_file, " in ",
        round(difftime(Sys.time(), start, units = "mins"), 1), " min")

# Arm-level values from the consensus GISTIC results (WGS only), used for the
# 1p/19q codeletion and monosomy 6 (6p and 6q) rules. Only the rows needed are kept.
arm_zip <- file.path(root_dir, "data", "cnv-consensus-gistic.zip")
arm_tmp <- tempfile("gistic")
dir.create(arm_tmp)
unzip(arm_zip, files = "cnv-consensus-gistic/broad_values_by_arm.txt",
      exdir = arm_tmp, junkpaths = TRUE)

arm_df <- data.table::fread(file.path(arm_tmp, "broad_values_by_arm.txt"),
                            data.table = FALSE) %>%
  filter(`Chromosome Arm` %in% c("1p", "1q", "19q", "6p", "6q", "22q", "17p", "17q", "12p", "9q")) %>%
  column_to_rownames("Chromosome Arm") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Kids_First_Biospecimen_ID")

# 1p and 19q are ploidy-aware: the GISTIC broad value is not corrected for the
# ploidy of the tumor (7316-4220, ploidy 3: 1p reads +1.33). For the WGS tumors with
# a ploidy estimate, the 1p and 19q values are the mean deviation of the consensus
# segment copy number from the tumor ploidy over the whole arm (copies; neutral
# segments count as ploidy; loss is <= -0.5, as for the GISTIC values); the other tumors keep the GISTIC
# value. The GISTIC values stay in 1p_gistic and 19q_gistic.
arm_regions <- tibble(
  arm = c("1p", "19q"), chrom = c("chr1", "chr19"),
  arm_start = c(0, 28700000), arm_end = c(123400000, 58617616)
)
tumor_ploidy <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                         guess_max = 100000, show_col_types = FALSE) %>%
  filter(experimental_strategy == "WGS", sample_type == "Tumor") %>%
  transmute(Kids_First_Biospecimen_ID, ploidy = suppressWarnings(as.numeric(tumor_ploidy))) %>%
  filter(!is.na(ploidy)) %>%
  distinct()
arm_cn <- read_tsv(file.path(root_dir, "data", "cnv-consensus.seg.gz"),
                   show_col_types = FALSE, col_types = cols(.default = "c")) %>%
  transmute(Kids_First_Biospecimen_ID = ID, chrom,
            loc.start = as.numeric(loc.start), loc.end = as.numeric(loc.end),
            copy.num = suppressWarnings(as.numeric(copy.num))) %>%
  filter(chrom %in% arm_regions$chrom, !is.na(copy.num)) %>%
  inner_join(arm_regions, by = "chrom", relationship = "many-to-many") %>%
  mutate(overlap = pmin(loc.end, arm_end) - pmax(loc.start, arm_start)) %>%
  filter(overlap > 0) %>%
  inner_join(tumor_ploidy, by = "Kids_First_Biospecimen_ID") %>%
  # neutral segments and gaps have no copy number in the consensus file, so they
  # count as the tumor ploidy: the deviation from ploidy is summed over the
  # called segments and divided by the whole arm length (a few focal losses
  # do not read as an arm loss)
  group_by(Kids_First_Biospecimen_ID, arm) %>%
  summarise(relative = sum((copy.num - ploidy) * overlap) /
              (first(arm_end) - first(arm_start)),
            .groups = "drop") %>%
  pivot_wider(names_from = arm, values_from = relative)
arm_df <- arm_df %>%
  rename(`1p_gistic` = `1p`, `19q_gistic` = `19q`) %>%
  left_join(arm_cn, by = "Kids_First_Biospecimen_ID") %>%
  mutate(`1p` = coalesce(`1p`, `1p_gistic`), `19q` = coalesce(`19q`, `19q_gistic`))

write_tsv(arm_df, file.path(scratch_dir, "arm_values.tsv"))
message("Wrote arm values for ", nrow(arm_df), " biospecimens")
unlink(arm_tmp, recursive = TRUE)

# Expression of the genes in input/expression_annotations.tsv (e.g. EZHIP
# overexpression), subset from the collapsed TPM matrix
expr_genes <- read_tsv(file.path(module_dir, "input", "expression_annotations.tsv"),
                       show_col_types = FALSE)
# genes of the annotations and of the expression (expr) lesion rules
expr_gene_names <- unique(c(expr_genes$gene, unlist(str_split(
  rules$gene[rules$lesion_type == "expr"], "\\|"))))
# genes are rows of the collapsed matrix, so the parquet copy answers this with
# a row filter; without it the RDS is read
tpm_parquet <- file.path(root_dir, "data", "parquet",
                         "gene-expression-rsem-tpm-collapsed.parquet")
if (file.exists(tpm_parquet)) {
  message("Reading gene-expression-rsem-tpm-collapsed from parquet")
  con_tpm <- connect_duckdb()
  tpm <- dbGetQuery(con_tpm, sprintf(
    "SELECT * FROM read_parquet('%s') WHERE gene_symbol IN (%s)",
    normalizePath(tpm_parquet),
    paste0("'", expr_gene_names, "'", collapse = ", "))) %>%
    tibble::column_to_rownames("gene_symbol")
  dbDisconnect(con_tpm, shutdown = TRUE)
} else {
  message("Reading gene-expression-rsem-tpm-collapsed.rds (no parquet copy)")
  tpm <- readRDS(file.path(root_dir, "data", "gene-expression-rsem-tpm-collapsed.rds"))
}
stopifnot(all(expr_gene_names %in% rownames(tpm)))

# Overexpression is a z-score of log2(TPM + 1) within RNA library type, against
# PBTA tumor RNA-seq biospecimens. Library type matters: exome-capture
# libraries have much higher EZHIP TPM than poly-A or stranded ones.
rna_meta <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                     guess_max = 100000, show_col_types = FALSE) %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  transmute(Kids_First_Biospecimen_ID, cohort, sample_type, pathology_diagnosis,
            library_group = case_when(
              str_detect(RNA_library, "exome") ~ "exome capture",
              str_detect(RNA_library, "poly-A") ~ "poly-A",
              TRUE ~ "stranded"
            ))

expr_long <- tpm[expr_gene_names, , drop = FALSE] %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "Kids_First_Biospecimen_ID", values_to = "tpm") %>%
  inner_join(rna_meta, by = "Kids_First_Biospecimen_ID") %>%
  mutate(log2_tpm = log2(tpm + 1))

# TBXT (brachyury) is expressed by chordomas, which would inflate the reference
# and hide the overexpression it is meant to show, so they are left out of its
# reference
reference <- expr_long %>%
  filter(cohort == "PBTA", sample_type == "Tumor",
         !(gene == "TBXT" & pathology_diagnosis %in% "Chordoma")) %>%
  group_by(gene, library_group) %>%
  summarise(ref_mean = mean(log2_tpm), ref_sd = sd(log2_tpm), n_ref = n(),
            .groups = "drop")

expr_long <- expr_long %>%
  left_join(reference, by = c("gene", "library_group")) %>%
  mutate(z = (log2_tpm - ref_mean) / ref_sd) %>%
  select(gene, Kids_First_Biospecimen_ID, library_group, tpm, z)

write_tsv(expr_long, file.path(scratch_dir, "expression_annotation_genes.tsv"))
message("Wrote expression for ", length(expr_gene_names), " gene(s) in ",
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
  cnv_out <- file.path(scratch_dir, "cnv_lesion_genes.parquet")
  cnv_query <- sprintf(
    "COPY (
       SELECT biospecimen_id AS Kids_First_Biospecimen_ID, gene_symbol, status
       FROM %s
       WHERE gene_symbol IN (%s) AND status NOT IN ('neutral', 'unchanged')
     ) TO '%s' (FORMAT PARQUET)",
    release_source("consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only",
                   "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz",
                   root_dir),
    paste0("'", cnv_genes, "'", collapse = ", "), cnv_out
  )
  con2 <- connect_duckdb()
  dbExecute(con2, cnv_query)
  dbDisconnect(con2, shutdown = TRUE)
  message("Wrote gene-level copy number for ", length(cnv_genes), " gene(s)")
}

# Structural variants (Manta) for the genes of sv rules (e.g. the BCOR exon 15
# internal tandem duplication, which the SNV callers miss)
sv_genes <- rules %>%
  filter(lesion_type == "sv") %>%
  pull(gene) %>%
  str_split("\\|") %>%
  unlist() %>%
  unique()

sv_out <- file.path(scratch_dir, "sv_lesion_genes.parquet")
con4 <- connect_duckdb()
dbExecute(con4, sprintf(
  "COPY (
     SELECT DISTINCT \"Kids.First.Biospecimen.ID.Tumor\" AS Kids_First_Biospecimen_ID,
            Gene_name, ALT || ' ' || coalesce(Location, '') AS sv_call,
            CAST(SV_length AS VARCHAR) AS sv_length
     FROM %s
     WHERE regexp_matches(Gene_name, '(^|;)(%s)(;|$)')
   ) TO '%s' (FORMAT PARQUET)",
  release_source("sv-manta", "sv-manta.tsv.gz", root_dir),
  paste(sv_genes, collapse = "|"), sv_out))
dbDisconnect(con4, shutdown = TRUE)
message("Wrote Manta structural variants for ", length(sv_genes), " gene(s)")

# LGG module lesions that need more than a gene list ------------------------------
# CDKN2A / CDKN2B copy number (legacy LGG rule: copy number <= 1)
cdkn_out <- file.path(scratch_dir, "cdkn2_copy_number.parquet")
con3 <- connect_duckdb()
dbExecute(con3, sprintf(
  "COPY (
     SELECT biospecimen_id AS Kids_First_Biospecimen_ID, gene_symbol,
            TRY_CAST(copy_number AS DOUBLE) AS copy_number, status
     FROM %s
     WHERE gene_symbol IN ('CDKN2A', 'CDKN2B')
   ) TO '%s' (FORMAT PARQUET)",
  release_source("consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only",
                 "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz",
                 root_dir),
  cdkn_out))

# FGFR1 kinase-domain duplication, as in the LGG module: an intragenic tandem
# duplication (Manta, 100 kb or shorter; the kinase domain duplication is about
# 5 kb, longer duplications are broad gains of the gene) that spans the FGFR1 protein tyrosine kinase domain, or a consensus
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
            SV_start, SV_end, ALT, TRY_CAST(SV_length AS DOUBLE) AS SV_length
     FROM %s
     WHERE Gene_name = 'FGFR1' AND ALT = '<DUP:TANDEM>'
       AND TRY_CAST(SV_length AS DOUBLE) <= 100000
   ) TO '%s' (FORMAT PARQUET)",
  release_source("sv-manta", "sv-manta.tsv.gz", root_dir), manta_out))

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

# LGG / HGG kinase fusions (reciprocal fusions included) and Ewing sarcoma
# fusions of interest, from the annoFuse putative oncogenic and DGD fusions
fusion_calls <- fusion_calls_of_interest(root_dir)
write_tsv(fusion_calls$lgg_hgg, file.path(scratch_dir, "fusion_calls_lgg_hgg.tsv"))
write_tsv(fusion_calls$ewings, file.path(scratch_dir, "fusion_calls_ewings.tsv"))
message("Wrote the fusion calls of interest for ", nrow(fusion_calls$lgg_hgg),
        " biospecimens")
