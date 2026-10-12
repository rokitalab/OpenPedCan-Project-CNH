################################################################################
# 05-compile-subtypes.R
#
# One row per match_id with the module's subtype and its evidence.
#
# Subtype precedence:
#   1. high-confidence methylation (03-validate-subtypes.R)
#   2. lesion-only fallback (03-validate-subtypes.R), when assigned
#   3. LGG lesion-level call, for tumors with no methylation or fallback call
#
# `lesion_refinement` adds the LGG lesion-level subtype (KIAA1549-BRAF, BRAF
# V600E, NF1, ...) on top of the methylation family when LGG made a specific
# call. TP53 columns are annotations and are not part of the subtype.
#
# Usage: Rscript --vanilla 05-compile-subtypes.R
################################################################################

suppressPackageStartupMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS")
results_dir <- file.path(module_dir, "results")

# free-text evidence columns are often entirely empty, so force them to character
evidence_types <- cols(tp53_snv_evidence = col_character(),
                       tp53_cnv_evidence = col_character(),
                       egfr_variants = col_character(),
                       .default = col_guess())

# chordoma: SMARCB1 focal copy status and expression, as annotations
smarcb1 <- read_tsv(
  file.path(root_dir, "analyses", "molecular-subtyping-chordoma", "results",
            "chordoma_smarcb1_status.tsv"),
  col_types = cols(.default = col_character())
) %>%
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(
    smarcb1_sample_id = first(sample_id),
    smarcb1_participant = first(Kids_First_Participant_ID),
    smarcb1_focal_status = first(na.omit(focal_SMARCB1_status)),
    smarcb1_expression = first(na.omit(SMARCB1_expression)),
    .groups = "drop"
  )

# manual review (input/manual_subtype_review.tsv): a subtype set by a pathology
# report or clinical review. It has the highest precedence; the automated call is
# kept in automated_subtype and disagreements with it are flagged.
manual_review <- read_tsv(file.path(module_dir, "input", "manual_subtype_review.tsv"),
                          show_col_types = FALSE) %>%
  transmute(sample_id, manual_subtype = molecular_subtype,
            manual_source = source, manual_note = note)
stopifnot(!anyDuplicated(manual_review$sample_id))

# umbrella groups for subtypes (input/subtype_groups.tsv)
subtype_groups <- read_tsv(file.path(module_dir, "input", "subtype_groups.tsv"),
                           show_col_types = FALSE)

umbrella_group <- function(subtype) {
  if (is.na(subtype)) return(NA_character_)
  hit <- subtype_groups %>% filter(str_detect(subtype, subtype_regex)) %>% slice(1)
  if (nrow(hit) == 0) NA_character_ else hit$subtype_group
}

# LGG labels that do not name a specific lesion
nonspecific <- "To be classified|wildtype|NOS"

lgg <- read_tsv(
  file.path(root_dir, "analyses", "molecular-subtyping-LGAT", "results",
            "lgat_subtyping.tsv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(match_id), !is.na(molecular_subtype)) %>%
  transmute(match_id, lgg_subtype = molecular_subtype) %>%
  distinct() %>%
  group_by(match_id) %>%
  summarise(lgg_subtype = paste(sort(unique(lgg_subtype)), collapse = " | "),
            .groups = "drop") %>%
  mutate(lgg_specific = !str_detect(lgg_subtype, nonspecific))

# medulloblastoma: expression classifier, concordance with methylation, and
# the SHH subtype (04-validate-mb.R)
mb <- read_tsv(file.path(results_dir, "mb_validation.tsv"),
               guess_max = 10000, show_col_types = FALSE) %>%
  filter(!is.na(match_id)) %>%
  select(match_id, mb_sample_id = sample_id, mb_expression_group,
         mb_concordance, mb_shh_subtype, mb_shh_source)

# one methylation call per tumor: highest score wins when several biospecimens
methyl <- read_tsv(file.path(results_dir, "methylation_subtype_validation.tsv"),
                   guess_max = 10000, col_types = evidence_types) %>%
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  arrange(desc(methyl_score), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(subtype_source = if_else(methylation_overridden %in% TRUE,
                                  "lesion overrides methylation", "methylation"))

fallback <- read_tsv(file.path(results_dir, "lesion_fallback_subtypes.tsv"),
                     guess_max = 10000, col_types = evidence_types) %>%
  filter(assign, !match_id %in% methyl$match_id) %>%
  mutate(subtype_source = "lesion fallback")

free_text <- read_tsv(file.path(results_dir, "free_text_dmg.tsv"),
                      guess_max = 10000, col_types = evidence_types) %>%
  filter(!match_id %in% c(methyl$match_id, fallback$match_id)) %>%
  mutate(subtype_source = "pathology free text")

lesion_negative <- read_tsv(file.path(results_dir, "lesion_negative_hgg.tsv"),
                            guess_max = 10000, col_types = evidence_types) %>%
  filter(!match_id %in% c(methyl$match_id, fallback$match_id)) %>%
  mutate(subtype_source = "lesion negative")

shared <- c("match_id", "Kids_First_Participant_ID", "sample_id",
            "OPC_molecular_subtype", "subtype_source", "evidence_level",
            "lesion_evidence", "pathology_diagnosis", "dx_group",
            "dx_concordance")

compiled <- bind_rows(
  methyl %>% select(all_of(shared), methyl_class, methyl_score, lesion_status,
                    lesion_contradiction,
                    methyl_conflict, starts_with("tp53_"),
                    ends_with("_mutation_status"), ends_with("_variants"),
                    ends_with("_amplification"), ends_with("_activating_alteration"),
                    ends_with("_tpm"), ends_with("_z"), ends_with("_expression_status"),
                    starts_with("h3_"), starts_with("h3k27_")),
  free_text %>% select(all_of(shared), starts_with("tp53_"),
                       ends_with("_mutation_status"), ends_with("_variants"),
                       ends_with("_amplification"),
                       ends_with("_activating_alteration"),
                       starts_with("h3_"), starts_with("h3k27_")),
  lesion_negative %>% select(all_of(shared), starts_with("tp53_"),
                             ends_with("_mutation_status"), ends_with("_variants"),
                             ends_with("_amplification"),
                             ends_with("_activating_alteration"),
                             starts_with("h3_")),
  fallback %>% select(all_of(shared), starts_with("tp53_"),
                      ends_with("_mutation_status"), ends_with("_variants"),
                    ends_with("_amplification"), ends_with("_activating_alteration"),
                    ends_with("_tpm"), ends_with("_z"), ends_with("_expression_status"),
                    starts_with("h3_"), starts_with("h3k27_"))
) %>%
  full_join(lgg, by = "match_id") %>%
  full_join(mb, by = "match_id") %>%
  # full join: chordoma tumors with no other source can still be subtyped
  # from focal SMARCB1 loss
  full_join(smarcb1, by = "match_id") %>%
  mutate(
    sample_id = coalesce(sample_id, mb_sample_id, smarcb1_sample_id),
    Kids_First_Participant_ID = coalesce(Kids_First_Participant_ID,
                                         smarcb1_participant),
    lesion_refinement = if_else(lgg_specific %in% TRUE, lgg_subtype,
                                NA_character_),
    # LGG is only the subtype when nothing else assigned one
    # expression classifier group, when it names a single group
    mb_expression_single = if_else(str_detect(replace_na(mb_expression_group, ""), "\\|"),
                                   NA_character_, mb_expression_group),
    subtype_source = case_when(
      !is.na(subtype_source) ~ subtype_source,
      lgg_specific %in% TRUE ~ "LGG lesion",
      !is.na(mb_expression_single) ~ "expression classifier",
      # chordoma diagnosis (the chordoma module's cohort) with focal SMARCB1 loss
      smarcb1_focal_status %in% "loss" ~ "SMARCB1 focal loss",
      TRUE ~ NA_character_
    ),
    molecular_subtype_cns = case_when(
      !is.na(OPC_molecular_subtype) ~ OPC_molecular_subtype,
      subtype_source == "LGG lesion" ~ lgg_subtype,
      subtype_source == "expression classifier" ~ mb_expression_single,
      subtype_source == "SMARCB1 focal loss" ~ "CHORDM, poorly differentiated"
    ),
    evidence_level = case_when(
      subtype_source == "pathology free text" ~ "free text",
      subtype_source %in% c("LGG lesion", "SMARCB1 focal loss") ~ "lesion only",
      subtype_source == "expression classifier" ~ "expression only",
      TRUE ~ evidence_level
    ),
    methyl_conflict = replace_na(methyl_conflict, FALSE)
  ) %>%
  # manual review overrides the automated subtype
  left_join(manual_review, by = "sample_id") %>%
  mutate(
    automated_subtype = molecular_subtype_cns,
    manual_conflict = !is.na(manual_subtype) & !is.na(automated_subtype) &
      manual_subtype != automated_subtype,
    subtype_source = if_else(!is.na(manual_subtype), "manual review", subtype_source),
    evidence_level = if_else(!is.na(manual_subtype), "manual review", evidence_level),
    molecular_subtype_cns = coalesce(manual_subtype, molecular_subtype_cns)
  ) %>%
  # umbrella group (e.g. the pediatric H3-wildtype HGG classes), unless an H3
  # variant was found, which contradicts the wildtype grouping
  mutate(subtype_group = map_chr(molecular_subtype_cns, umbrella_group)) %>%
  mutate(
    h3_group_conflict = !is.na(subtype_group) & h3_status %in% "mutant",
    subtype_group = case_when(
      h3_group_conflict ~ molecular_subtype_cns,
      !is.na(subtype_group) ~ subtype_group,
      TRUE ~ molecular_subtype_cns
    ),
    needs_review = is.na(manual_subtype) &
      (dx_concordance %in% "discordant" | methyl_conflict |
         h3_group_conflict | h3_variant_conflict %in% TRUE |
         mb_concordance %in% "discordant" | !is.na(lesion_contradiction))
  ) %>%
  select(match_id, Kids_First_Participant_ID, sample_id,
         molecular_subtype_cns, subtype_group, lesion_refinement, subtype_source,
         automated_subtype, manual_conflict, manual_source, manual_note,
         evidence_level, lesion_status, lesion_evidence, lesion_contradiction,
         methyl_class,
         mb_expression_group, mb_concordance, mb_shh_subtype,
         smarcb1_focal_status, smarcb1_expression,
         methyl_score, methyl_conflict, pathology_diagnosis, dx_group,
         dx_concordance, needs_review, starts_with("tp53_"),
         ends_with("_mutation_status"), ends_with("_variants"),
         ends_with("_amplification"), ends_with("_activating_alteration"),
         ends_with("_tpm"), ends_with("_z"), ends_with("_expression_status"),
         starts_with("h3_"), starts_with("h3k27_")) %>%
  arrange(match_id)

stopifnot(!anyDuplicated(compiled$match_id))
write_tsv(compiled, file.path(results_dir, "molecular_subtype_cns.tsv"))

message("Tumors with a subtype: ",
        sum(!is.na(compiled$molecular_subtype_cns)), " of ", nrow(compiled))
print(count(compiled, subtype_source))
message("With lesion refinement: ", sum(!is.na(compiled$lesion_refinement)))
