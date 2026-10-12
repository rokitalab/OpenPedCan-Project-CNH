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
scratch_dir <- file.path(root_dir, "scratch", "molecular-subtyping-CNS")
dir.create(scratch_dir, showWarnings = FALSE, recursive = TRUE)

maf_file <- file.path(root_dir, "data", "snv-consensus-plus-hotspots.maf.tsv.gz")
out_file <- file.path(scratch_dir, "snv_lesion_genes.parquet")

rules <- read_tsv(file.path(module_dir, "input", "subtype_lesion_rules.tsv"),
                  show_col_types = FALSE)

annotation_genes <- read_tsv(file.path(module_dir, "input", "annotation_genes.tsv"),
                             show_col_types = FALSE)

genes <- rules %>%
  filter(lesion_type == "snv") %>%
  pull(gene) %>%
  str_split("\\|") %>%
  unlist() %>%
  c(annotation_genes$gene) %>%
  unique()

gene_sql <- paste0("'", genes, "'", collapse = ", ")

con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

# all_varchar and a short sample avoid type sniffing over the whole file
query <- sprintf(
  "COPY (
     SELECT Hugo_Symbol, HGVSp_Short, Variant_Classification,
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
  con2 <- dbConnect(duckdb::duckdb())
  dbExecute(con2, cnv_query)
  dbDisconnect(con2, shutdown = TRUE)
  message("Wrote gene-level copy number for ", length(cnv_genes), " gene(s)")
}
