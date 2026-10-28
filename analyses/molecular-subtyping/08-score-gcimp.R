################################################################################
# 08-score-gcimp.R
#
# Glioma CpG island methylator phenotype (G-CIMP) score of every methylation
# array: the mean beta value over the G-CIMP probes that are hypermethylated in
# IDH-mutant astrocytoma. The 1,503 probes are the G-CIMP probe set of
# Noushmehr et al. 2010 (Cancer Cell 17:510-522, Figure 1; the supplementary
# file on https://gdc.cancer.gov/about-data/publications/gbm_2010), in
# input/gcimp_probes.tsv; the published call is a clustering of those probes in
# the TCGA arrays, which is not reproduced here.
#
# The direction of each probe is learned from this cohort: a probe is used when
# its mean beta is more than 0.2 higher in the IDH-mutant astrocytoma arrays
# (DKFZ A_IDH_LG, A_IDH_HG) than in the arrays of non-IDH gliomas (pediatric
# high-grade glioma, DMG, pilocytic astrocytoma, PXA and GBM classes). The score
# is evidence for an IDH-mutant glioma methylation class, not a call: the
# percentile column is the fraction of the non-IDH glioma arrays with a lower
# score. The reference arrays are scored with the rest, so their scores are
# slightly optimistic.
#
# Needs the beta value matrix (data/parquet/methyl-beta-values.parquet; bash
# scripts/download-parquet.sh methyl-beta-values). Run after 07-compile-subtypes.R
# (it reads the methylation calls) and before 09-compare-to-legacy.R.
#
# Usage: Rscript --vanilla 08-score-gcimp.R
# Output: results/gcimp_scores.tsv
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
source(file.path(module_dir, "utils", "duckdb-connect.R"))

mut_classes <- c("A_IDH_LG", "A_IDH_HG")
nonidh_classes <- c("DMG_K27", "DMG_EGFR", "PA_INF", "PA_MID", "PA_CORT", "PXA",
                    "pedHGG_A", "pedHGG_RTK1A", "pedHGG_RTK1B", "pedHGG_RTK1C",
                    "pedHGG_RTK2A", "pedHGG_RTK2B", "pedHGG_MYCN", "HGG_E", "HGAP",
                    "GBM_MES_TYP", "GBM_MES_ATYP", "GBM_RTK1")
min_difference <- 0.2

probes <- read_tsv(file.path(module_dir, "input", "gcimp_probes.tsv"),
                   show_col_types = FALSE)$ProbeID

all_calls <- read_tsv(file.path(module_dir, "results", "methylation_calls.tsv"),
                      show_col_types = FALSE) %>%
  distinct(Kids_First_Biospecimen_ID, methyl_source, methyl_class, methyl_score)
# the reference arrays are the DKFZ v12.8 calls; the class shown is the call used
# for the subtype (DKFZ, or NIH when DKFZ is below its threshold)
calls <- all_calls %>% filter(methyl_source == "dkfz_v12.8")

### Beta values of the G-CIMP probes -------------------------------------------

con <- connect_duckdb()
beta_file <- file.path(root_dir, "data", "parquet", "methyl-beta-values.parquet")
stopifnot(file.exists(beta_file))
arrays <- dbGetQuery(con, sprintf(
  "SELECT column_name FROM (DESCRIBE SELECT * FROM read_parquet('%s'))",
  beta_file))$column_name[-1]
beta <- dbGetQuery(con, sprintf(
  "SELECT Probe_ID, %s FROM read_parquet('%s') WHERE Probe_ID IN (%s)",
  paste(sprintf('"%s"', arrays), collapse = ", "), beta_file,
  paste(sprintf("'%s'", probes), collapse = ", ")))
dbDisconnect(con, shutdown = TRUE)

# some probes are listed more than once: their rows are averaged
mat <- do.call(rbind, lapply(split(seq_len(nrow(beta)), beta$Probe_ID), function(i) {
  colMeans(as.matrix(beta[i, -1, drop = FALSE]), na.rm = TRUE)
}))
message("G-CIMP probes on the array: ", nrow(mat), " of ", length(probes))

### Probes hypermethylated in IDH-mutant astrocytoma ---------------------------

ids_of <- function(classes) {
  calls %>% filter(methyl_class %in% classes, Kids_First_Biospecimen_ID %in% arrays) %>%
    pull(Kids_First_Biospecimen_ID)
}
mut_ids <- ids_of(mut_classes)
nonidh_ids <- ids_of(nonidh_classes)
difference <- rowMeans(mat[, mut_ids, drop = FALSE], na.rm = TRUE) -
  rowMeans(mat[, nonidh_ids, drop = FALSE], na.rm = TRUE)
hyper <- names(difference)[which(difference > min_difference)]
message("Reference arrays: ", length(mut_ids), " IDH-mutant astrocytoma, ",
        length(nonidh_ids), " non-IDH glioma; hypermethylated probes: ", length(hyper))

score <- colMeans(mat[hyper, , drop = FALSE], na.rm = TRUE)
nonidh_scores <- score[nonidh_ids]

hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  filter(Kids_First_Biospecimen_ID %in% arrays) %>%
  distinct(Kids_First_Biospecimen_ID, sample_id, match_id)

out <- tibble(Kids_First_Biospecimen_ID = names(score), gcimp_score = unname(score),
              gcimp_probes = length(hyper)) %>%
  mutate(gcimp_percentile = map_dbl(gcimp_score, ~ mean(nonidh_scores < .x)),
         gcimp_reference = case_when(
           Kids_First_Biospecimen_ID %in% mut_ids ~ "IDH-mutant astrocytoma",
           Kids_First_Biospecimen_ID %in% nonidh_ids ~ "non-IDH glioma",
           TRUE ~ NA_character_)) %>%
  left_join(all_calls %>% filter(!is.na(methyl_source)) %>%
              group_by(Kids_First_Biospecimen_ID) %>%
              slice(1) %>% ungroup(),
            by = "Kids_First_Biospecimen_ID") %>%
  left_join(hist, by = "Kids_First_Biospecimen_ID") %>%
  select(Kids_First_Biospecimen_ID, sample_id, match_id, methyl_source, methyl_class, methyl_score,
         gcimp_score, gcimp_percentile, gcimp_probes, gcimp_reference) %>%
  arrange(desc(gcimp_score))

write_tsv(out, file.path(module_dir, "results", "gcimp_scores.tsv"))
print(out %>% group_by(gcimp_reference) %>%
        summarise(n = n(), median = median(gcimp_score), min = min(gcimp_score),
                  max = max(gcimp_score)))
