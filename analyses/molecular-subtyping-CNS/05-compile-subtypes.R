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
                       EGFR_variants = col_character(),
                       .default = col_guess())

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

# one methylation call per tumor: highest score wins when several biospecimens
methyl <- read_tsv(file.path(results_dir, "methylation_subtype_validation.tsv"),
                   guess_max = 10000, col_types = evidence_types) %>%
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  arrange(desc(methyl_score), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(subtype_source = "methylation")

fallback <- read_tsv(file.path(results_dir, "lesion_fallback_subtypes.tsv"),
                     guess_max = 10000, col_types = evidence_types) %>%
  filter(assign, !match_id %in% methyl$match_id) %>%
  mutate(subtype_source = "lesion fallback")

shared <- c("match_id", "Kids_First_Participant_ID", "sample_id",
            "OPC_molecular_subtype", "subtype_source", "evidence_level",
            "lesion_evidence", "pathology_diagnosis", "dx_group",
            "dx_concordance")

compiled <- bind_rows(
  methyl %>% select(all_of(shared), methyl_class, methyl_score, lesion_status,
                    methyl_conflict, starts_with("tp53_"),
                    ends_with("_mutation_status"), ends_with("_variants"),
                    ends_with("_amplification"), ends_with("_activating_alteration"),
                    starts_with("h3k27_")),
  fallback %>% select(all_of(shared), starts_with("tp53_"),
                      ends_with("_mutation_status"), ends_with("_variants"),
                    ends_with("_amplification"), ends_with("_activating_alteration"),
                    starts_with("h3k27_"))
) %>%
  full_join(lgg, by = "match_id") %>%
  mutate(
    lesion_refinement = if_else(lgg_specific %in% TRUE, lgg_subtype,
                                NA_character_),
    # LGG is only the subtype when nothing else assigned one
    subtype_source = case_when(
      !is.na(subtype_source) ~ subtype_source,
      lgg_specific %in% TRUE ~ "LGG lesion",
      TRUE ~ NA_character_
    ),
    molecular_subtype_cns = case_when(
      !is.na(OPC_molecular_subtype) ~ OPC_molecular_subtype,
      subtype_source == "LGG lesion" ~ lgg_subtype
    ),
    evidence_level = if_else(subtype_source == "LGG lesion", "lesion only",
                             evidence_level),
    methyl_conflict = replace_na(methyl_conflict, FALSE),
    needs_review = dx_concordance %in% "discordant" | methyl_conflict
  ) %>%
  select(match_id, Kids_First_Participant_ID, sample_id,
         molecular_subtype_cns, lesion_refinement, subtype_source,
         evidence_level, lesion_status, lesion_evidence, methyl_class,
         methyl_score, methyl_conflict, pathology_diagnosis, dx_group,
         dx_concordance, needs_review, starts_with("tp53_"),
         ends_with("_mutation_status"), ends_with("_variants"),
                    ends_with("_amplification"), ends_with("_activating_alteration"),
                    starts_with("h3k27_")) %>%
  arrange(match_id)

stopifnot(!anyDuplicated(compiled$match_id))
write_tsv(compiled, file.path(results_dir, "molecular_subtype_cns.tsv"))

message("Tumors with a subtype: ",
        sum(!is.na(compiled$molecular_subtype_cns)), " of ", nrow(compiled))
print(count(compiled, subtype_source))
message("With lesion refinement: ", sum(!is.na(compiled$lesion_refinement)))
