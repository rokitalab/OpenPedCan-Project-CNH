################################################################################
# 06-subtype-solid-tumors.R
#
# Subtypes of the solid tumors in the cohort: neuroblastic tumors
# (neuroblastoma, ganglioneuroblastoma, ganglioneuroma) by MYCN status.
# Ported from molecular-subtyping-NBL.
#
# MYCN evidence, in this order:
#   1. copy number: the consensus gene-level call "amplification" on a DNA
#      assay of the tumor decides (reviewer decision)
#   2. the clinically reported MYCN status (the free-text diagnosis, filled from
#      the TARGET and GMKF clinical MYCN tables when it is missing), which
#      outranks expression (reviewer decision, 7316 NLF cell line)
#   3. expression: MYCN TPM at or above the cutoff (140.83, the legacy module's
#      cutoff from the TPM distribution of the NBL samples) on an RNA assay;
#      it also decides when the copy number says not amplified and the report
#      says amplified
#   4. a copy number call other than amplification
# A tumor whose call differs from the clinical report is flagged. A tumor with a
# copy number that says not amplified and a clinical report that says
# amplified, and no expression, is left without a subtype and flagged.
#
# Ewing sarcoma (ported from molecular-subtyping-EWS): the label is Rare CNS, EWS
# (CNS cohorts; EWS for the extracranial PPTC xenografts) for a
# tumor with a high-confidence EWS methylation class, a hallmark EWSR1 / FUS
# fusion in the RNA (EWSR1--FLI1, EWSR1--ERG, FUS--ERG, from the fusion
# summary) in a CNS cohort, or an Ewing sarcoma diagnosis. Methylation comes
# first, then the fusion, then the diagnosis alone; a methylation class or a
# fusion in a tumor with another diagnosis, and a diagnosis with a different
# methylation class, are flagged.
#
# Output: results/solid_subtypes.tsv (one row per match_id)
#
# Usage: Rscript --vanilla 06-subtype-solid-tumors.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
input_dir <- file.path(module_dir, "input")
results_dir <- file.path(module_dir, "results")
source(file.path(module_dir, "utils", "duckdb-connect.R"))
source(file.path(module_dir, "utils", "release-tables.R"))
source(file.path(module_dir, "utils", "dx-concordance.R"))
source(file.path(module_dir, "utils", "mismatch-arrays.R"))

# MYCN TPM cutoff from the legacy molecular-subtyping-NBL module
mycn_tpm_cutoff <- 140.83

### Tumors -----------------------------------------------------------------------

nbl_dx <- jsonlite::fromJSON(
  file.path(input_dir, "path_dx_strings.json"))$nbl$exact_path_dx

hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  filter(sample_type == "Tumor", !is.na(match_id))

nbl_hist <- hist %>% filter(pathology_diagnosis %in% nbl_dx)

# clinical MYCN status, to fill free-text diagnoses that do not state it
clinical <- bind_rows(
  read_tsv(file.path(input_dir, "nbl_target_clinical_mycn_status.tsv"),
           show_col_types = FALSE) %>%
    transmute(Kids_First_Participant_ID = USI_ID, clinical = TARGET_MYCN_STATUS),
  read_tsv(file.path(input_dir, "nbl_gmkf_clinical_mycn_status.tsv"),
           show_col_types = FALSE) %>%
    transmute(Kids_First_Participant_ID = USI_ID, clinical = GMKF_MYCN_STATUS)
) %>%
  distinct(Kids_First_Participant_ID, .keep_all = TRUE)

# clinical MYCN status of the PPTC neuroblastoma models, keyed on the model name
# (sample_id); optional: input/nbl_pptc_clinical_mycn_status.tsv with the columns
# sample_id and PPTC_MYCN_STATUS ("MYCN amp" or "MYCN non-amp")
pptc_mycn_file <- file.path(input_dir, "nbl_pptc_clinical_mycn_status.tsv")
clinical_by_sample <- if (file.exists(pptc_mycn_file)) {
  read_tsv(pptc_mycn_file, show_col_types = FALSE) %>%
    transmute(sample_id, clinical_sample = PPTC_MYCN_STATUS) %>%
    filter(!is.na(clinical_sample)) %>%
    distinct(sample_id, .keep_all = TRUE)
} else {
  tibble(sample_id = character(), clinical_sample = character())
}

status_from_text <- function(x) {
  case_when(
    str_detect(x, regex("non.?amp|not amplified", ignore_case = TRUE)) ~ "non-amplified",
    str_detect(x, regex("mycn[^|]*amp", ignore_case = TRUE)) ~ "amplified",
    TRUE ~ NA_character_
  )
}

tumors <- nbl_hist %>%
  group_by(match_id) %>%
  summarise(
    Kids_First_Participant_ID = first(Kids_First_Participant_ID),
    sample_id = first(sample_id),
    cohort = first(cohort),
    composition = first(composition),
    pathology_diagnosis = first(pathology_diagnosis),
    free_text = paste(unique(na.omit(pathology_free_text_diagnosis)), collapse = " | "),
    .groups = "drop"
  ) %>%
  left_join(clinical, by = "Kids_First_Participant_ID") %>%
  left_join(clinical_by_sample, by = "sample_id") %>%
  mutate(
    clinical = coalesce(clinical, clinical_sample),
    mycn_text_status = coalesce(status_from_text(free_text), status_from_text(clinical)),
    mycn_text_source = case_when(
      !is.na(status_from_text(free_text)) ~ "free-text diagnosis",
      !is.na(status_from_text(clinical)) ~ "clinical MYCN table",
      TRUE ~ NA_character_
    )
  )

dna_ids <- nbl_hist %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Sequencing"),
         is.na(RNA_library)) %>%
  distinct(match_id, id = Kids_First_Biospecimen_ID)
rna_ids <- nbl_hist %>%
  filter(experimental_strategy == "RNA-Seq" |
           (experimental_strategy == "Targeted Sequencing" & !is.na(RNA_library))) %>%
  distinct(match_id, id = Kids_First_Biospecimen_ID)

### MYCN copy number and expression ----------------------------------------------

con <- connect_duckdb()
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

id_list <- function(ids) paste0("'", ids, "'", collapse = ", ")

cn <- dbGetQuery(con, sprintf(
  "SELECT biospecimen_id AS id, status FROM %s
   WHERE gene_symbol = 'MYCN' AND biospecimen_id IN (%s)",
  release_source("consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only",
                 "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz",
                 root_dir),
  id_list(dna_ids$id))) %>%
  as_tibble() %>%
  inner_join(dna_ids, by = "id") %>%
  mutate(status = str_to_lower(status)) %>%
  group_by(match_id) %>%
  summarise(mycn_cn_status = if_else(any(status == "amplification"),
                                     "amplification",
                                     paste(sort(unique(status)), collapse = ";")),
            .groups = "drop")

tpm_parquet <- file.path(root_dir, "data", "parquet",
                         "gene-expression-rsem-tpm-collapsed.parquet")
if (file.exists(tpm_parquet)) {
  message("Reading gene-expression-rsem-tpm-collapsed from parquet")
  tpm_wide <- dbGetQuery(con, sprintf(
    "SELECT * FROM read_parquet('%s') WHERE gene_symbol = 'MYCN'",
    normalizePath(tpm_parquet))) %>%
    as_tibble()
} else {
  message("Reading gene-expression-rsem-tpm-collapsed.rds (no parquet copy)")
  m <- readRDS(file.path(root_dir, "data", "gene-expression-rsem-tpm-collapsed.rds"))
  tpm_wide <- as_tibble(m["MYCN", , drop = FALSE], rownames = "gene_symbol")
}
mycn_tpm <- tpm_wide %>%
  pivot_longer(-gene_symbol, names_to = "id", values_to = "tpm") %>%
  inner_join(rna_ids, by = "id") %>%
  group_by(match_id) %>%
  summarise(mycn_tpm = max(tpm, na.rm = TRUE), .groups = "drop")

### Methylation (annotation, and a conflict check) -------------------------------

methyl <- read_tsv(file.path(results_dir, "methylation_calls.tsv"),
                   show_col_types = FALSE) %>%
  exclude_mismatched_arrays(results_dir) %>%
  filter(methyl_confidence == "high") %>%
  inner_join(hist %>% filter(experimental_strategy == "Methylation") %>%
               distinct(Kids_First_Biospecimen_ID, match_id),
             by = "Kids_First_Biospecimen_ID") %>%
  group_by(match_id) %>%
  summarise(methyl_class = paste(unique(methyl_class), collapse = ";"),
            .groups = "drop")

### Subtype ----------------------------------------------------------------------

solid <- tumors %>%
  left_join(cn, by = "match_id") %>%
  left_join(mycn_tpm, by = "match_id") %>%
  left_join(methyl, by = "match_id") %>%
  mutate(
    tpm_call = case_when(is.na(mycn_tpm) ~ NA_character_,
                         mycn_tpm >= mycn_tpm_cutoff ~ "amplified",
                         TRUE ~ "non-amplified"),
    cn_call = case_when(is.na(mycn_cn_status) ~ NA_character_,
                        mycn_cn_status == "amplification" ~ "amplified",
                        TRUE ~ "non-amplified"),
    # copy number amplification decides; otherwise the clinical report, except
    # when the copy number says not amplified and the report says amplified,
    # which expression decides when it exists; then expression; then a copy
    # number call that is not amplification
    mycn_call = case_when(
      cn_call %in% "amplified" ~ "amplified",
      !is.na(mycn_text_status) &
        !(cn_call %in% "non-amplified" & mycn_text_status == "amplified") ~
        mycn_text_status,
      !is.na(tpm_call) ~ tpm_call,
      cn_call %in% "non-amplified" & mycn_text_status %in% "amplified" ~ NA_character_,
      !is.na(cn_call) ~ "non-amplified",
      TRUE ~ NA_character_
    ),
    # the evidence that supports the call
    subtype_source = case_when(
      is.na(mycn_call) ~ NA_character_,
      !is.na(cn_call) & cn_call == mycn_call ~ "MYCN copy number",
      !is.na(mycn_text_status) & mycn_text_status == mycn_call ~ "clinical MYCN status",
      !is.na(tpm_call) & tpm_call == mycn_call ~ "MYCN expression",
      TRUE ~ NA_character_
    ),
    molecular_subtype = case_when(
      mycn_call == "amplified" ~ "NBL, MYCN amplified",
      mycn_call == "non-amplified" ~ "NBL, MYCN non-amplified",
      TRUE ~ NA_character_
    ),
    evidence_level = case_when(
      is.na(molecular_subtype) ~ NA_character_,
      subtype_source == "clinical MYCN status" ~ "clinical report",
      subtype_source == "MYCN expression" ~ "expression only",
      TRUE ~ "lesion only"
    ),
    subtype_alterations_confirmed = case_when(
      mycn_call == "amplified" & subtype_source == "MYCN copy number" ~
        "MYCN amplification",
      mycn_call == "amplified" & subtype_source == "MYCN expression" ~
        "MYCN overexpression",
      TRUE ~ NA_character_
    ),
    # the clinical report disagrees with the molecular call, or the evidence
    # conflicts and no call was made
    text_conflict = !is.na(mycn_call) & !is.na(mycn_text_status) &
      mycn_call != mycn_text_status & subtype_source != "clinical MYCN status",
    evidence_conflict = is.na(mycn_call) & cn_call %in% "non-amplified" &
      mycn_text_status %in% "amplified",
    methyl_conflict = !is.na(methyl_class) & !is.na(molecular_subtype) &
      ((str_detect(methyl_class, "NB_MYCN") & mycn_call == "non-amplified")),
    # a high-confidence methylation class of another tumor type (a mix-up of
    # the array, or a misdiagnosis); the MYCN call is kept
    methyl_other_class = !is.na(methyl_class) &
      !str_detect(methyl_class, "^(NB_|CNS_NB|CTRL_|INFLAM_)"),
    subtype_review_flag = case_when(
      evidence_conflict ~ "MYCN copy number not amplified, clinical report amplified",
      text_conflict ~ "MYCN call differs from the clinical report",
      methyl_conflict ~ "NB_MYCN methylation, MYCN not amplified",
      methyl_other_class ~ "methylation class is not a neuroblastic class",
      TRUE ~ NA_character_
    ),
    Notes = case_when(
      subtype_source == "MYCN copy number" ~ "Subtype from MYCN copy number",
      subtype_source == "MYCN expression" ~ "Subtype from MYCN expression",
      subtype_source == "clinical MYCN status" ~ "Subtype from clinical MYCN status",
      TRUE ~ NA_character_
    )
  ) %>%
  select(match_id, Kids_First_Participant_ID, sample_id, cohort, composition,
         pathology_diagnosis, molecular_subtype, subtype_source, evidence_level,
         subtype_alterations_confirmed, subtype_review_flag, Notes,
         mycn_cn_status, mycn_tpm, mycn_text_status, mycn_text_source, methyl_class)

### Ewing sarcoma ----------------------------------------------------------------

ews_dx <- jsonlite::fromJSON(
  file.path(input_dir, "path_dx_strings.json"))$ews$exact_path_dx

methyl_all <- read_tsv(file.path(results_dir, "methylation_calls.tsv"),
                       show_col_types = FALSE) %>%
  exclude_mismatched_arrays(results_dir) %>%
  filter(methyl_confidence == "high") %>%
  inner_join(hist %>% filter(experimental_strategy == "Methylation") %>%
               distinct(Kids_First_Biospecimen_ID, match_id),
             by = "Kids_First_Biospecimen_ID") %>%
  group_by(match_id) %>%
  summarise(methyl_class = paste(unique(methyl_class), collapse = ";"),
            .groups = "drop")

ews_foi <- read_tsv(file.path(root_dir, "scratch", "molecular-subtyping", "fusion_calls_ewings.tsv"),
                    show_col_types = FALSE)
fusion_cols <- setdiff(names(ews_foi), "Kids_First_Biospecimen_ID")
ews_fusions <- ews_foi %>%
  pivot_longer(all_of(fusion_cols), names_to = "fusion", values_to = "n") %>%
  filter(n > 0) %>%
  inner_join(hist %>% distinct(Kids_First_Biospecimen_ID, match_id, cohort),
             by = "Kids_First_Biospecimen_ID") %>%
  # the CNS cohort (PBTA, which includes its DGD sub-cohort) and the PPTC
  # xenografts; FUS--ERG is also an AML fusion
  filter(cohort %in% c("PBTA", "PPTC")) %>%
  group_by(match_id) %>%
  summarise(ews_fusions = paste(sort(unique(fusion)), collapse = "; "),
            .groups = "drop")

ews_tumors <- hist %>%
  group_by(match_id) %>%
  summarise(
    Kids_First_Participant_ID = first(Kids_First_Participant_ID),
    sample_id = first(sample_id),
    cohort = first(cohort),
    composition = first(composition),
    pathology_diagnosis = first(pathology_diagnosis),
    .groups = "drop"
  ) %>%
  left_join(methyl_all, by = "match_id") %>%
  left_join(ews_fusions, by = "match_id") %>%
  mutate(is_ews_dx = pathology_diagnosis %in% ews_dx,
         methyl_ews = str_detect(replace_na(methyl_class, ""), "(^|;)EWS($|;)")) %>%
  filter(is_ews_dx | methyl_ews | !is.na(ews_fusions))

ews <- ews_tumors %>%
  mutate(
    # a high-confidence methylation class of another tumor type outranks the
    # diagnosis alone (the class sets the subtype and the Ewing diagnosis is
    # flagged); a hallmark fusion still gives EWS
    other_class_dx_only = is_ews_dx & !methyl_ews & is.na(ews_fusions) &
      !is.na(methyl_class) & !str_detect(methyl_class, "^(CTRL_|INFLAM_)"),
    # Ewing sarcoma of the CNS cohorts is a rare CNS tumor; the extracranial
    # xenograft models (PPTC) are EWS
    molecular_subtype = case_when(
      other_class_dx_only ~ NA_character_,
      # PBTA includes the DGD and Kentucky sub-cohorts (CNS); the DGD cohort is
      # not CNS
      cohort %in% "PBTA" ~ "Rare CNS, EWS",
      TRUE ~ "EWS"
    ),
    subtype_source = case_when(
      other_class_dx_only ~ NA_character_,
      methyl_ews ~ "methylation",
      !is.na(ews_fusions) ~ "EWSR1 fusion",
      TRUE ~ "pathology diagnosis"
    ),
    evidence_level = case_when(
      other_class_dx_only ~ NA_character_,
      methyl_ews ~ "methylation",
      !is.na(ews_fusions) ~ "lesion only",
      TRUE ~ "diagnosis"
    ),
    subtype_alterations_confirmed = ews_fusions,
    # a methylation class of another tumor type in an Ewing diagnosis (a
    # mismatched array, e.g. 7316-1856), or EWS evidence in another diagnosis
    other_class = is_ews_dx & !methyl_ews & !is.na(methyl_class) &
      !str_detect(methyl_class, "^(CTRL_|INFLAM_)"),
    other_dx = !is_ews_dx & (methyl_ews | !is.na(ews_fusions)),
    subtype_review_flag = case_when(
      other_class ~ "methylation class is not EWS",
      other_dx ~ "EWS evidence in another diagnosis",
      TRUE ~ NA_character_
    ),
    Notes = case_when(
      is.na(molecular_subtype) ~ NA_character_,
      subtype_source == "methylation" ~ "Subtype from methylation",
      subtype_source == "EWSR1 fusion" ~ "Subtype from EWSR1 fusion",
      TRUE ~ "Subtype from pathology diagnosis"
    )
  ) %>%
  select(match_id, Kids_First_Participant_ID, sample_id, cohort, composition,
         pathology_diagnosis, molecular_subtype, subtype_source, evidence_level,
         subtype_alterations_confirmed, subtype_review_flag, Notes,
         methyl_class, ews_fusions)

# the neuroblastic diagnosis wins for a tumor in both sets (none expected)
ews <- ews %>% filter(!match_id %in% solid$match_id)
solid <- bind_rows(solid, ews)

# concordance with the pathology diagnosis, as for the other subtypes
dx_groups <- read_tsv(file.path(input_dir, "pathology_dx_groups.tsv"),
                      show_col_types = FALSE)
family_compat <- read_tsv(file.path(input_dir, "subtype_family_dx_compatibility.tsv"),
                          show_col_types = FALSE) %>%
  arrange(order)
solid <- solid %>%
  mutate(dx_group = dx_to_groups(pathology_diagnosis, dx_groups),
         dx_concordance = dx_concordance(molecular_subtype, dx_group, family_compat))

stopifnot(!anyDuplicated(solid$match_id))
write_tsv(solid, file.path(results_dir, "solid_subtypes.tsv"))

message("Solid tumors: ", nrow(solid), " (", nrow(ews), " Ewing sarcoma); with a subtype: ",
        sum(!is.na(solid$molecular_subtype)))
print(count(solid, molecular_subtype, subtype_source))
print(count(solid, subtype_review_flag))
