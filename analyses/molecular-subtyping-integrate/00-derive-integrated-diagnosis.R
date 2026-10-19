################################################################################
# 00-derive-integrated-diagnosis.R
#
# Derive integrated_diagnosis, broad_histology and short_histology for the
# subtypes of molecular-subtyping by rule (WHO CNS 2021 naming, keeping the
# legacy strings where they exist):
#   1. an exact label in input/subtype_diagnosis_overrides.tsv, otherwise
#   2. the family default in input/subtype_diagnosis_families.tsv:
#      "<integrated_prefix>, <the label after its first comma>".
# The first comma-separated token of integrated_diagnosis is the cancer_group
# in 01-integrate-subtyping.Rmd, so family prefixes are the existing cancer
# group names.
#
# Writes results/integrated_diagnosis.tsv (one row per tumor) and
# results/subtype_diagnosis_map.tsv (one row per label, with the legacy value
# for comparison).
#
# Usage: Rscript 00-derive-integrated-diagnosis.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-integrate")
cns_file <- file.path(root_dir, "analyses", "molecular-subtyping", "results",
                      "molecular_subtype.tsv")
legacy_file <- file.path(root_dir, "analyses", "archive", "molecular-subtyping-pathology",
                         "results",
                         "compiled_molecular_subtypes_with_clinical_pathology_feedback_and_report_info.tsv")
results_dir <- file.path(module_dir, "results")
dir.create(results_dir, showWarnings = FALSE)

families <- read_tsv(file.path(module_dir, "input", "subtype_diagnosis_families.tsv"),
                     show_col_types = FALSE, na = character()) %>%
  mutate(across(c(broad_histology, short_histology), ~ na_if(.x, "")))
overrides <- read_tsv(file.path(module_dir, "input", "subtype_diagnosis_overrides.tsv"),
                      show_col_types = FALSE, na = character()) %>%
  mutate(across(c(broad_histology, short_histology), ~ na_if(.x, "")))
stopifnot(!anyDuplicated(overrides$molecular_subtype),
          !anyDuplicated(families$family))

cns <- read_tsv(cns_file, show_col_types = FALSE, guess_max = 100000)

labels <- cns %>%
  filter(!is.na(molecular_subtype)) %>%
  count(molecular_subtype, name = "tumors") %>%
  mutate(family = str_extract(molecular_subtype, "^[^,]*"),
         rest = str_trim(str_remove(molecular_subtype, "^[^,]*,?"))) %>%
  left_join(families, by = "family") %>%
  left_join(overrides, by = "molecular_subtype", suffix = c("", "_override")) %>%
  mutate(
    integrated_diagnosis = case_when(
      !is.na(integrated_diagnosis) ~ integrated_diagnosis,
      !is.na(integrated_prefix) & rest != "" ~ paste0(integrated_prefix, ", ", rest),
      !is.na(integrated_prefix) ~ integrated_prefix
    ),
    rule = case_when(
      molecular_subtype %in% overrides$molecular_subtype ~ "override",
      !is.na(integrated_prefix) ~ "family",
      TRUE ~ "none"
    ),
    broad_histology = coalesce(broad_histology_override, broad_histology),
    short_histology = coalesce(short_histology_override, short_histology)
  ) %>%
  select(molecular_subtype, tumors, rule, integrated_diagnosis,
         broad_histology, short_histology)

# compare with the legacy pathology table (most frequent pair per label)
legacy <- read_tsv(legacy_file, show_col_types = FALSE, guess_max = 100000) %>%
  filter(!is.na(molecular_subtype), molecular_subtype != "NA") %>%
  count(molecular_subtype, integrated_diagnosis, broad_histology, short_histology) %>%
  group_by(molecular_subtype) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(molecular_subtype = molecular_subtype,
         legacy_integrated_diagnosis = integrated_diagnosis,
         legacy_broad_histology = broad_histology,
         legacy_short_histology = short_histology)

map <- labels %>%
  left_join(legacy, by = "molecular_subtype") %>%
  mutate(differs_from_legacy = !is.na(legacy_integrated_diagnosis) &
           integrated_diagnosis != legacy_integrated_diagnosis) %>%
  arrange(molecular_subtype)
write_tsv(map, file.path(results_dir, "subtype_diagnosis_map.tsv"))

# per tumor: the label's diagnosis, except that the neuroblastic tumors take
# the diagnosis name of the tumor (neuroblastoma, ganglioneuroblastoma,
# ganglioneuroma) with the MYCN status, as in the legacy module
per_tumor <- cns %>%
  filter(!is.na(molecular_subtype)) %>%
  select(match_id, Kids_First_Participant_ID, sample_id, molecular_subtype,
         subtype_source, subtype_review_flag, subtype_alterations_confirmed,
         alteration_confirmed, Notes, methyl_class, pathology_diagnosis) %>%
  left_join(select(map, molecular_subtype, integrated_diagnosis,
                   broad_histology, short_histology),
            by = "molecular_subtype") %>%
  mutate(
    nbl_name = str_trim(str_extract(replace_na(pathology_diagnosis, ""), "^[^,]*")),
    is_nbl = str_detect(molecular_subtype, "^NBL, MYCN") &
      str_detect(nbl_name, regex("neuroblast|ganglioneuro", ignore_case = TRUE)),
    integrated_diagnosis = if_else(
      is_nbl,
      paste0(nbl_name, ", ",
             if_else(str_detect(molecular_subtype, "non-amplified"),
                     "MYCN non-amplified", "MYCN-amplified")),
      integrated_diagnosis
    ),
    short_histology = if_else(is_nbl, nbl_name, short_histology)
  ) %>%
  select(-nbl_name, -is_nbl, -pathology_diagnosis)
write_tsv(per_tumor, file.path(results_dir, "integrated_diagnosis.tsv"))

message("Labels: ", nrow(map), "; no rule: ", sum(map$rule == "none"),
        "; differ from legacy: ", sum(map$differs_from_legacy))
print(map %>% filter(rule == "none" | differs_from_legacy) %>%
        select(molecular_subtype, tumors, rule, integrated_diagnosis,
               legacy_integrated_diagnosis), n = 100)
