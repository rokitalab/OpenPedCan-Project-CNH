################################################################################
# tumor-evidence.R
#
# Everything the module holds about one or more tumors, in one place, for a
# review: the clinical record, the DKFZ v12.8 and NIH v3.1 methylation calls,
# the subtype row, the variants, fusions and copy number of the lesion genes,
# the 1p/19q arm values and the G-CIMP score. It reads the tables of
# run-molecular-subtyping.sh (results/ and scratch/molecular-subtyping/), so run
# the module first.
#
# Usage: Rscript --vanilla utils/tumor-evidence.R 7316-6389 7316-526
#   (sample_id or match_id; run from analyses/molecular-subtyping)
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
results_dir <- file.path(module_dir, "results")
scratch_dir <- file.path(root_dir, "scratch", "molecular-subtyping")
source(file.path(module_dir, "utils", "duckdb-connect.R"))

ids <- commandArgs(trailingOnly = TRUE)
if (length(ids) == 0) stop("Give one or more sample_id or match_id values")

show <- function(title, x) {
  cat("\n--", title, "--\n")
  if (is.null(x) || nrow(x) == 0) cat("(none)\n") else print(as.data.frame(x), row.names = FALSE, right = FALSE)
}

hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  filter(sample_id %in% ids | match_id %in% ids)
if (nrow(hist) == 0) stop("No tumor with these ids in histologies-base.tsv")
bs_ids <- hist$Kids_First_Biospecimen_ID

con <- connect_duckdb()
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
parquet <- function(name, where) {
  dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s') WHERE %s",
                          file.path(scratch_dir, name), where)) %>% as_tibble()
}
bs_in <- paste(sprintf("'%s'", bs_ids), collapse = ", ")

subtype <- read_tsv(file.path(results_dir, "molecular_subtype.tsv"),
                    guess_max = 100000, show_col_types = FALSE,
                    col_types = cols(.default = col_guess())) %>%
  filter(sample_id %in% ids | match_id %in% ids)
bti <- read_tsv(file.path(module_dir, "input", "bti_methylation_classifications.tsv"),
                show_col_types = FALSE, col_types = cols(.default = "c")) %>%
  filter(Bioassay_ID %in% bs_ids)

for (id in unique(hist$sample_id)) {
  cat("\n==============================", id, "==============================\n")
  h <- hist %>% filter(sample_id == id)
  ids_of <- h$Kids_First_Biospecimen_ID

  show("Clinical", h %>%
         distinct(match_id, cohort, tumor_descriptor, CNS_region, primary_site,
                  age_years = round(age_at_diagnosis_days / 365.25, 1),
                  pathology_diagnosis, pathology_free_text_diagnosis,
                  cancer_predispositions))
  show("Assays", h %>% count(experimental_strategy, RNA_library))

  show("Methylation (DKFZ v12.8 and NIH v3.1 subclass, with scores)", bti %>%
         filter(Bioassay_ID %in% ids_of) %>%
         transmute(Bioassay_ID,
                   dkfz_v12.8 = paste0(`dkfz_v12.8_subclass`, " ",
                                       signif(as.numeric(`dkfz_v12.8_subclass_maxscore`), 2)),
                   nih_v3.1 = paste0(NIH_v3.1_subclass, " ",
                                     signif(as.numeric(NIH_v3.1_subclass_meanscore), 2))))

  s <- subtype %>% filter(sample_id == id)
  show("Subtype row", s %>%
         select(molecular_subtype, subtype_source, lesion_refinement,
                subtype_alterations_confirmed, subtype_review_flag,
                methylation_mismatch, diagnosis_review, dx_concordance,
                lesion_status, lesion_evidence, lesion_contradiction, Notes) %>%
         mutate(across(everything(), as.character)) %>%
         pivot_longer(everything(), names_to = "field", values_to = "value") %>%
         filter(!is.na(value)))

  snv <- parquet("snv_lesion_genes.parquet",
                 sprintf("Kids_First_Biospecimen_ID IN (%s) AND HGVSp_Short IS NOT NULL",
                         bs_in)) %>%
    filter(Kids_First_Biospecimen_ID %in% ids_of, !str_detect(HGVSp_Short, "=$"))
  show("Variants in the lesion genes (protein changes)",
       snv %>% distinct(Hugo_Symbol, HGVSp_Short, Variant_Classification) %>%
         arrange(Hugo_Symbol))

  onc <- read_tsv(file.path(root_dir, "data", "fusion-putative-oncogenic.tsv"),
                  show_col_types = FALSE, col_select = c(Sample, FusionName, Fusion_Type)) %>%
    filter(Sample %in% ids_of)
  show("Putative oncogenic fusions", onc %>% distinct(FusionName, Fusion_Type))

  show("Copy number of the lesion genes",
       parquet("cnv_lesion_genes.parquet",
               sprintf("Kids_First_Biospecimen_ID IN (%s)", bs_in)) %>%
         filter(Kids_First_Biospecimen_ID %in% ids_of) %>%
         distinct(gene_symbol, status) %>% arrange(gene_symbol))

  arms <- read_tsv(file.path(scratch_dir, "arm_values.tsv"), show_col_types = FALSE) %>%
    filter(Kids_First_Biospecimen_ID %in% ids_of)
  show("Arm values (1p and 19q are ploidy-aware; loss <= -0.5)", arms)

  gc_file <- file.path(results_dir, "gcimp_scores.tsv")
  if (file.exists(gc_file)) {
    show("G-CIMP score (08-score-gcimp.R)",
         read_tsv(gc_file, show_col_types = FALSE) %>%
           filter(Kids_First_Biospecimen_ID %in% ids_of) %>%
           select(Kids_First_Biospecimen_ID, gcimp_score, gcimp_percentile))
  }
}
