################################################################################
# make-methyl-testing-parquet.R
#
# Write a small testing copy of the methylation beta values for the CI jobs of
# the molecular subtyping, from the release parquet (data/<release>/parquet/
# methyl-beta-values.parquet, about 15 GB, probes in rows and arrays in columns).
#
# The subtyping reads the beta values in two places only, so the testing file has
# only what they use:
#   10-plot-mb-methylation.R   the arrays of the medulloblastoma tumors (all of
#                              them: UMAP and PCA of the 20,000 most variable
#                              probes of these arrays)
#   08-score-gcimp.R           the G-CIMP probes of Noushmehr et al. 2010
#                              (input/gcimp_probes.tsv) in the arrays of
#                              IDH-mutant astrocytoma (all), a sample of each
#                              non-IDH glioma class, and the tumors reviewed
#                              for G-CIMP
# Rows: the G-CIMP probes and the 20,000 probes of highest variance among the
# medulloblastoma arrays, so the two scripts give the same result as on the full
# file for the medulloblastoma plots. Columns: the arrays above.
#
# Output: data/<testing release>/parquet/methyl-beta-values.parquet (the same
# name as the release file, so scripts/download-parquet.sh fetches it from the
# testing release), plus its line in manifest.tsv and md5sum-parquet.txt. Upload
# the folder to the testing release in the bucket.
#
# Usage: Rscript scripts/make-methyl-testing-parquet.R [source release] [testing release]
#   (defaults v15 and testing). Needs analyses/molecular-subtyping/results/
#   methylation_calls.tsv (run 01-build-methylation-calls.R first).
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
  library(duckdb)
})

args <- commandArgs(trailingOnly = TRUE)
source_release <- if (length(args) >= 1) args[1] else "v15"
testing_release <- if (length(args) >= 2) args[2] else "testing"
set.seed(2026)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
source_file <- file.path(root_dir, "data", source_release, "parquet", "methyl-beta-values.parquet")
out_dir <- file.path(root_dir, "data", testing_release, "parquet")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "methyl-beta-values.parquet")
stopifnot(file.exists(source_file))

# per-glioma-class sample, and the tumors reviewed for G-CIMP
per_class <- 6
reviewed_samples <- c("7316-6389", "7316-526", "7316-1946", "7316-804", "7316-4917",
                      "7316-4993", "7316-4215", "7316-3165", "7316-4310")
mut_classes <- c("A_IDH_LG", "A_IDH_HG")
nonidh_classes <- c("DMG_K27", "DMG_EGFR", "PA_INF", "PA_MID", "PA_CORT", "PXA",
                    "pedHGG_A", "pedHGG_RTK1A", "pedHGG_RTK1B", "pedHGG_RTK1C",
                    "pedHGG_RTK2A", "pedHGG_RTK2B", "pedHGG_MYCN", "HGG_E", "HGAP",
                    "GBM_MES_TYP", "GBM_MES_ATYP", "GBM_RTK1")
n_variable <- 20000

con <- dbConnect(duckdb::duckdb(config = list(memory_limit = "8GB",
                                              temp_directory = file.path(root_dir, "scratch"))))
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
arrays_in_file <- dbGetQuery(con, sprintf(
  "SELECT column_name FROM (DESCRIBE SELECT * FROM read_parquet('%s'))", source_file
))$column_name[-1]

hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                 col_select = c(Kids_First_Biospecimen_ID, sample_id, match_id,
                                experimental_strategy, pathology_diagnosis),
                 show_col_types = FALSE) %>%
  filter(experimental_strategy == "Methylation",
         Kids_First_Biospecimen_ID %in% arrays_in_file)
calls <- read_tsv(file.path(module_dir, "results", "methylation_calls.tsv"),
                  show_col_types = FALSE) %>%
  inner_join(hist, by = "Kids_First_Biospecimen_ID")

### Arrays -----------------------------------------------------------------------

mb_ids <- calls %>%
  filter(str_detect(replace_na(OPC_molecular_subtype, ""), "^MB, ") |
           str_detect(replace_na(pathology_diagnosis, ""), regex("medulloblastoma", ignore_case = TRUE))) %>%
  pull(Kids_First_Biospecimen_ID)
mb_ids <- union(mb_ids, hist$Kids_First_Biospecimen_ID[
  str_detect(replace_na(hist$pathology_diagnosis, ""), regex("medulloblastoma", ignore_case = TRUE))])

mut_ids <- calls %>% filter(methyl_class %in% mut_classes) %>% pull(Kids_First_Biospecimen_ID)
nonidh_ids <- calls %>% filter(methyl_class %in% nonidh_classes) %>%
  group_by(methyl_class) %>% slice_sample(n = per_class) %>% ungroup() %>%
  pull(Kids_First_Biospecimen_ID)
reviewed_ids <- hist %>% filter(sample_id %in% reviewed_samples) %>%
  pull(Kids_First_Biospecimen_ID)

arrays <- unique(c(mb_ids, mut_ids, nonidh_ids, reviewed_ids))
message("Arrays: ", length(arrays), " (medulloblastoma ", length(mb_ids),
        ", IDH-mutant astrocytoma ", length(mut_ids), ", non-IDH glioma sample ",
        length(nonidh_ids), ", reviewed ", length(reviewed_ids), ")")

### Probes: G-CIMP probes and the most variable probes of the medulloblastoma arrays

gcimp <- read_tsv(file.path(module_dir, "input", "gcimp_probes.tsv"),
                  show_col_types = FALSE)$ProbeID
message("Reading the medulloblastoma arrays for the probe variance")
mb_cols <- paste(sprintf('"%s"', mb_ids), collapse = ", ")
total_probes <- dbGetQuery(con, sprintf(
  "SELECT count(*) AS n FROM read_parquet('%s')", source_file))$n
chunk <- 100000
variances <- vector("list", ceiling(total_probes / chunk))
for (i in seq_along(variances)) {
  x <- dbGetQuery(con, sprintf(
    "SELECT Probe_ID, %s FROM read_parquet('%s') LIMIT %d OFFSET %d",
    mb_cols, source_file, chunk, (i - 1) * chunk))
  m <- as.matrix(x[, -1])
  variances[[i]] <- tibble(Probe_ID = x$Probe_ID,
                           variance = apply(m, 1, var, na.rm = TRUE))
}
variable_probes <- bind_rows(variances) %>%
  filter(!is.na(variance)) %>%
  arrange(desc(variance)) %>%
  slice_head(n = n_variable) %>%
  pull(Probe_ID)
probes <- union(variable_probes, gcimp)
message("Probes: ", length(probes), " (", length(variable_probes), " most variable, ",
        length(gcimp), " G-CIMP)")

### Write ------------------------------------------------------------------------

sel <- paste(sprintf('"%s"', arrays), collapse = ", ")
probe_list <- paste(sprintf("'%s'", probes), collapse = ", ")
dbExecute(con, sprintf(
  "COPY (SELECT Probe_ID, %s FROM read_parquet('%s') WHERE Probe_ID IN (%s) ORDER BY Probe_ID)
   TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 2000)",
  sel, source_file, probe_list, out_file))

md5 <- unname(tools::md5sum(out_file))
size <- file.size(out_file)
message("Wrote ", out_file, " (", round(size / 1e6, 1), " MB, md5 ", md5, ")")

# manifest.tsv and md5sum-parquet.txt of the testing release: this table only (the
# rows are the selected probes, so they do not match the source rows)
manifest_file <- file.path(out_dir, "manifest.tsv")
source_rows <- dbGetQuery(con, sprintf(
  "SELECT count(*) AS n FROM read_parquet('%s')", source_file))$n
manifest <- tibble(table = "methyl-beta-values", source = basename(source_file),
                   rows = length(probes), source_rows = source_rows,
                   columns = length(arrays) + 1, size_mb = round(size / 1e6, 1),
                   md5 = md5, row_count_matches = FALSE)
write_tsv(manifest, manifest_file)
writeLines(paste(md5, "methyl-beta-values.parquet"), file.path(out_dir, "md5sum-parquet.txt"))
