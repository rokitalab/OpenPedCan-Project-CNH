################################################################################
# 02-prepare-lesion-tables.R
#
# Subset the SNV consensus + hotspot MAF to the genes named in
# input/subtype_lesion_rules.tsv and write a small parquet file, so that
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

genes <- rules %>%
  filter(lesion_type == "snv") %>%
  pull(gene) %>%
  str_split("\\|") %>%
  unlist() %>%
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
