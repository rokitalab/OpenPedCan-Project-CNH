################################################################################
# 03-validate-subtypes.R
#
# Validate methylation-derived subtypes (results/methylation_calls.tsv) against
#   1. genetic lesions (input/subtype_lesion_rules.tsv)
#   2. the historical pathology diagnosis (input/pathology_dx_groups.tsv and
#      input/subtype_family_dx_compatibility.tsv)
#
# Methylation is never overridden: lesions raise confidence, and a diagnosis
# that does not fit the subtype is flagged as a probable reclassification.
# PBTA tumors with no high-confidence methylation call are subtyped from
# lesions alone for rules marked use_as_fallback.
#
# Usage: Rscript --vanilla 03-validate-subtypes.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS")
input_dir <- file.path(module_dir, "input")
results_dir <- file.path(module_dir, "results")
scratch_dir <- file.path(root_dir, "scratch", "molecular-subtyping-CNS")

source(file.path(module_dir, "utils", "lesion-hits.R"))
source(file.path(module_dir, "utils", "dx-concordance.R"))
source(file.path(module_dir, "utils", "tp53-annotation.R"))
source(file.path(module_dir, "utils", "gene-annotation.R"))

### Inputs ---------------------------------------------------------------------

calls <- read_tsv(file.path(results_dir, "methylation_calls.tsv"),
                  show_col_types = FALSE)
rules <- read_tsv(file.path(input_dir, "subtype_lesion_rules.tsv"),
                  show_col_types = FALSE) %>%
  mutate(rule_order = row_number(),
         use_as_fallback = as.logical(use_as_fallback))
dx_groups <- read_tsv(file.path(input_dir, "pathology_dx_groups.tsv"),
                      show_col_types = FALSE)
family_compat <- read_tsv(
  file.path(input_dir, "subtype_family_dx_compatibility.tsv"),
  show_col_types = FALSE
) %>%
  arrange(order)

hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID, sample_id,
         match_id, cohort, sample_type, experimental_strategy,
         pathology_diagnosis)

snv_file <- file.path(scratch_dir, "snv_lesion_genes.parquet")
if (!file.exists(snv_file)) {
  stop("Run 02-prepare-lesion-tables.R first: ", snv_file, " not found")
}
con <- dbConnect(duckdb::duckdb())
snv <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", snv_file)) %>%
  as_tibble() %>%
  filter(!is.na(HGVSp_Short)) %>%
  distinct()
dbDisconnect(con, shutdown = TRUE)

fusion <- read_tsv(file.path(root_dir, "data", "fusion-annoFuse.tsv.gz"),
                   show_col_types = FALSE) %>%
  distinct(Kids_First_Biospecimen_ID = Sample, FusionName)

arm_values <- read_tsv(file.path(scratch_dir, "arm_1p19q.tsv"),
                       show_col_types = FALSE)

tp53 <- read_tp53_annotation(root_dir)
annotation_genes <- read_tsv(file.path(input_dir, "annotation_genes.tsv"),
                             show_col_types = FALSE)

### Lesion hits ----------------------------------------------------------------
# methylation is keyed on its own biospecimen, while SNVs, arm values and
# fusions are keyed on the DNA / RNA biospecimens of the same tumor, so they
# are linked through match_id

assay_ids <- function(strategies) {
  hist %>%
    filter(experimental_strategy %in% strategies, !is.na(match_id)) %>%
    distinct(match_id, id = Kids_First_Biospecimen_ID)
}
dna_ids <- assay_ids(c("WGS", "WXS", "Targeted Sequencing"))
rna_ids <- assay_ids("RNA-Seq")

lesion_hits <- find_lesion_hits(
  rules, snv, fusion, arm_values,
  assay_ids = list(snv = dna_ids, arm = dna_ids, fusion = rna_ids)
)

### Methylation-based subtypes -------------------------------------------------

rule_assays <- rules %>%
  group_by(OPC_molecular_subtype) %>%
  summarise(needs_dna = any(lesion_type %in% c("snv", "arm")),
            needs_rna = any(lesion_type == "fusion"),
            has_rule = TRUE, .groups = "drop")

confirming <- lesion_hits %>%
  group_by(match_id, OPC_molecular_subtype) %>%
  summarise(lesion_evidence = paste(sort(unique(evidence)), collapse = "; "),
            .groups = "drop")

validated <- calls %>%
  filter(methyl_confidence == "high", !methyl_control,
         !is.na(OPC_molecular_subtype)) %>%
  left_join(hist %>% select(Kids_First_Biospecimen_ID,
                            Kids_First_Participant_ID, sample_id, match_id,
                            pathology_diagnosis),
            by = "Kids_First_Biospecimen_ID") %>%
  left_join(rule_assays, by = "OPC_molecular_subtype") %>%
  left_join(confirming, by = c("match_id", "OPC_molecular_subtype")) %>%
  mutate(
    assay_available = (needs_dna %in% TRUE & match_id %in% dna_ids$match_id) |
      (needs_rna %in% TRUE & match_id %in% rna_ids$match_id),
    lesion_status = case_when(
      !(has_rule %in% TRUE) ~ "no_rule",
      !is.na(lesion_evidence) ~ "confirmed",
      assay_available ~ "not_detected",
      TRUE ~ "no_assay"
    ),
    # PXA (issue 113): DKFZ subclass PXA at score >= 0.8 is split by BRAF V600E
    # (the PXA rule's lesion). Without a DNA assay the status is not assessable.
    is_pxa = OPC_molecular_subtype == "PXA" &
      methyl_source %in% c("dkfz_v12.8", "dkfz_v12") & methyl_score >= 0.8,
    OPC_molecular_subtype = case_when(
      is_pxa & lesion_status == "confirmed" ~ "LGG, PXA BRAF V600E",
      is_pxa & lesion_status == "not_detected" ~ "LGG, PXA BRAF V600E negative",
      # PXA from a lower-confidence or NIH-only call is not split
      is_pxa | OPC_molecular_subtype == "PXA" ~ "LGG, PXA",
      TRUE ~ OPC_molecular_subtype
    ),
    dx_group = dx_to_groups(pathology_diagnosis, dx_groups),
    dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                    family_compat),
    evidence_level = if_else(lesion_status == "confirmed",
                             "methylation + lesion", "methylation only"),
    needs_review = dx_concordance == "discordant" | methyl_conflict
  ) %>%
  select(-is_pxa) %>%
  left_join(tp53, by = "match_id") %>%
  annotate_genes(annotation_genes, snv, dna_ids) %>%
  select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID, sample_id,
         match_id, OPC_molecular_subtype, methyl_class, methyl_source,
         methyl_score, methyl_conflict, methyl_conflicting_sources,
         lesion_status, lesion_evidence, pathology_diagnosis, dx_group,
         dx_concordance, starts_with("tp53_"),
         ends_with("_mutation_status"), ends_with("_variants"),
         evidence_level, needs_review)

write_tsv(validated, file.path(results_dir, "methylation_subtype_validation.tsv"))
validated %>%
  filter(needs_review) %>%
  arrange(dx_concordance, OPC_molecular_subtype) %>%
  write_tsv(file.path(results_dir, "subtypes_for_review.tsv"))

### Lesion-only fallback -------------------------------------------------------
# PBTA only (IDH1 R132 in TARGET leukemia is not an HGG), and not when the
# methylation call is a control / low-purity tissue, whose lesions are unreliable

control_match_ids <- calls %>%
  filter(methyl_confidence == "high", methyl_control) %>%
  inner_join(select(hist, Kids_First_Biospecimen_ID, match_id),
             by = "Kids_First_Biospecimen_ID") %>%
  pull(match_id) %>%
  unique()
pbta_match_ids <- unique(hist$match_id[hist$cohort == "PBTA"])

fallback <- lesion_hits %>%
  filter(rule_id %in% rules$rule_id[rules$use_as_fallback],
         match_id %in% pbta_match_ids,
         !match_id %in% validated$match_id,
         !match_id %in% control_match_ids) %>%
  left_join(select(rules, rule_id, rule_order), by = "rule_id") %>%
  group_by(match_id) %>%
  arrange(rule_order, evidence, .by_group = TRUE) %>%
  summarise(OPC_molecular_subtype = first(OPC_molecular_subtype),
            lesion_evidence = paste(unique(evidence), collapse = "; "),
            .groups = "drop") %>%
  # a match_id can carry more than one diagnosis (e.g. cell-line rows), which
  # are joined with ";" as in the histologies
  left_join(hist %>%
              filter(sample_type == "Tumor") %>%
              group_by(match_id) %>%
              summarise(Kids_First_Participant_ID = first(Kids_First_Participant_ID),
                        sample_id = first(sample_id),
                        pathology_diagnosis = paste(sort(unique(na.omit(pathology_diagnosis))),
                                                    collapse = ";"),
                        .groups = "drop") %>%
              mutate(pathology_diagnosis = na_if(pathology_diagnosis, "")),
            by = "match_id") %>%
  mutate(dx_group = dx_to_groups(pathology_diagnosis, dx_groups))

# IDH-mutant gliomas: with 1p/19q codeletion the tumor is an oligodendroglioma
# (Oligo, IDH). Otherwise the grade comes from the diagnosis group, and a tumor
# with an oligodendroglioma diagnosis but no codeletion is held unclassified.
codeleted <- unique(lesion_hits$match_id[lesion_hits$rule_id == "OLIGO_1P19Q"])
arm_assessed <- unique(dna_ids$match_id[dna_ids$id %in%
                                          arm_values$Kids_First_Biospecimen_ID])

fallback <- fallback %>%
  mutate(
    is_idh = OPC_molecular_subtype == "HGG, IDH",
    is_codeleted = match_id %in% codeleted,
    oligo_dx = str_detect(replace_na(pathology_diagnosis, ""),
                          regex("oligodendroglioma", ignore_case = TRUE)),
    hold_reason = case_when(
      is_idh & !is_codeleted & oligo_dx & match_id %in% arm_assessed ~
        "oligodendroglioma diagnosis, 1p/19q not codeleted",
      is_idh & !is_codeleted & oligo_dx ~
        "oligodendroglioma diagnosis, 1p/19q not assessable (no WGS arm data)"
    ),
    OPC_molecular_subtype = case_when(
      is_idh & is_codeleted ~ "Oligo, IDH",
      is_idh & !oligo_dx & dx_group == "LGG" ~ "LGG, IDH",
      TRUE ~ OPC_molecular_subtype
    ),
    lesion_evidence = if_else(is_idh & is_codeleted,
                              paste0(lesion_evidence, "; 1p+19q codeleted"),
                              lesion_evidence),
    dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                    family_compat),
    evidence_level = "lesion only",
    # a call the diagnosis contradicts, or a held call, is not assigned
    assign = dx_concordance != "discordant" & is.na(hold_reason)
  ) %>%
  select(-is_idh, -is_codeleted, -oligo_dx) %>%
  left_join(tp53, by = "match_id") %>%
  annotate_genes(annotation_genes, snv, dna_ids)

write_tsv(fallback, file.path(results_dir, "lesion_fallback_subtypes.tsv"))

### Summary --------------------------------------------------------------------

message("\n== Methylation-based subtypes: ", nrow(validated))
print(count(validated, lesion_status))
print(count(validated, dx_concordance))
message("Flagged for review: ", sum(validated$needs_review))

message("\n== Lesion-only fallback: ", nrow(fallback), " tumors without methylation")
print(count(fallback, OPC_molecular_subtype, assign, hold_reason))

# gaps in the lookup tables
gaps <- bind_rows(validated, fallback) %>%
  filter(dx_concordance == "subtype_family_unmapped" |
           str_detect(replace_na(dx_group, ""), "UNMAPPED"))
if (nrow(gaps) > 0) {
  message("\nSubtypes or diagnoses missing from the dx lookup tables:")
  print(distinct(gaps, OPC_molecular_subtype, pathology_diagnosis))
}
