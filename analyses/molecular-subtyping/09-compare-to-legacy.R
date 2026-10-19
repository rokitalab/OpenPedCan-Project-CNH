################################################################################
# 09-compare-to-legacy.R
#
# Compare molecular_subtype.tsv (07-compile-subtypes.R) with the subtypes of the
# previous release (the molecular_subtype column of data/histologies.tsv, made by
# the legacy per-diagnosis modules), per match_id. The legacy ", TP53" suffix is dropped, since TP53 is now an
# annotation, and DMG labels with `-mutant`, `-altered`, `EGFR` or `EZHIP` terms and `DHG, H3 G35-mutant` are compared as the legacy `DMG, H3 K28` and `DHG, H3 G35`.
# A legacy lesion-level label that equals `lesion_refinement` counts
# as agreement with the methylation family.
#
# Usage: Rscript --vanilla 09-compare-to-legacy.R
################################################################################

suppressPackageStartupMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
results_dir <- file.path(root_dir, "analyses", "molecular-subtyping",
                         "results")

# legacy labels renamed in this module to carry an HGG / LGG prefix
legacy_renames <- c("IHG, NTRK-altered" = "HGG, IHG, NTRK-altered",
                    "IHG, ROS1-altered" = "HGG, IHG, ROS1-altered",
                    "IHG, ALK-altered" = "HGG, IHG, ALK-altered",
                    "IHG, MET-altered" = "HGG, IHG, MET-altered",
                    "IHG, To be classified" = "HGG, IHG",
                    "HGAP" = "HGG, HGAP", "HGNET" = "HGG, HGNET",
                    "IHG" = "HGG, IHG", "SEGA" = "LGG, SEGA",
                    "GBM, MES" = "HGG, MES", "GBM, RTK" = "HGG, RTK")

# legacy labels that mean "no specific subtype was assigned"
# "ETMR, NOS" is a defined subtype (lesion and diagnosis based), not a default
unassigned <- "To be classified|(?<!ETMR), NOS$|wildtype$"

legacy <- read_tsv(
  file.path(root_dir, "data", "histologies.tsv"),
  col_select = c(match_id, sample_type, molecular_subtype),
  col_types = cols(.default = col_character())
) %>%
  filter(sample_type == "Tumor", !is.na(match_id), !is.na(molecular_subtype)) %>%
  # one label per tumor: the biospecimens of a match_id share it
  summarise(molecular_subtype = paste(sort(unique(molecular_subtype)), collapse = " | "),
            .by = match_id) %>%
  transmute(match_id,
            legacy_subtype = str_remove(molecular_subtype, ", TP53$") %>%
              # labels renamed here to carry an HGG / LGG / GNT prefix
              recode(!!!legacy_renames),
            legacy_group = str_extract(legacy_subtype, "^[^,]+")) %>%
  distinct()

new <- read_tsv(file.path(results_dir, "molecular_subtype.tsv"),
                guess_max = 10000, col_types = cols(.default = col_guess(),
                  tp53_snv_evidence = col_character(),
                  tp53_cnv_evidence = col_character(),
                  egfr_variants = col_character())) %>%
  # the legacy DMG, H3 K28 and DHG, H3 G35 cover both the methylation-only
  # label and the "-mutant" label (variant confirmed) used here
  mutate(new_subtype = molecular_subtype %>%
           str_replace("^DMG, H3 K28(-mutant|-altered)?(, EGFR)?(, EZHIP)?$",
                       "DMG, H3 K28") %>%
           str_remove("(?<=H3 G35)-mutant$")) %>%
  select(match_id, sample_id, new_subtype, subtype_source, lesion_refinement)

comparison <- full_join(legacy, new, by = "match_id") %>%
  mutate(
    legacy_assigned = !is.na(legacy_subtype) &
      !str_detect(legacy_subtype, unassigned),
    outcome = case_when(
      is.na(legacy_subtype) & is.na(new_subtype) ~ "neither",
      is.na(legacy_subtype) ~ "new only",
      is.na(new_subtype) & legacy_assigned ~ "legacy only (specific)",
      is.na(new_subtype) ~ "legacy only (unassigned)",
      legacy_subtype == new_subtype ~ "same",
      # the CNS Ewing sarcoma label of the release is "EWS" (vocabulary)
      legacy_subtype == "EWS" & new_subtype == "Rare CNS, EWS" ~ "same (vocabulary)",
      legacy_subtype %in% lesion_refinement ~ "same (lesion refinement)",
      # the legacy label is the new one without its lesion: the legacy base (a
      # trailing ", NOS" removed) followed by the lesion, e.g. ETMR, NOS ->
      # ETMR, C19MC-altered or PXA -> PXA BRAF V600E
      str_starts(new_subtype, fixed(paste0(str_remove(legacy_subtype, ", NOS$"), ", "))) |
        str_starts(new_subtype, fixed(paste0(str_remove(legacy_subtype, ", NOS$"), " "))) ~
        "same (lesion refinement)",
      !legacy_assigned ~ "newly assigned",
      # the Rare CNS subtypes were added to the repository after the release: not
      # compared
      str_starts(new_subtype, "Rare CNS") | str_starts(legacy_subtype, "Rare CNS") ~
        "Rare CNS (added after the release)",
      TRUE ~ "different"
    )
  )

# reviewed column: JLR = reviewed by the reviewer (initials), only for the rows
# whose outcome is "different" (and the Rare CNS rows). Reviewed as a group: a MYCN subtype from the
# clinical report, a MYCN subtype from expression or copy number (the release
# predates the NBL MYCN TPM fixes of 2026-10-06 and 2026-10-08) and a manual review
# (a reviewer decision); the other tumors are
# listed in input/legacy_comparison_review.tsv with a note.
reviewed_tumors <- read_tsv(file.path(root_dir, "analyses", "molecular-subtyping", "input",
                                      "legacy_comparison_review.tsv"),
                            show_col_types = FALSE) %>%
  distinct(sample_id, .keep_all = TRUE)
comparison <- comparison %>%
  left_join(reviewed_tumors %>% select(sample_id, reviewed_listed = reviewed,
                                       note_listed = note), by = "sample_id") %>%
  mutate(
    # only the tumors whose outcome is "different" are reviewed, and the Rare CNS
    # subtypes (subtyped and reviewed by the reviewer when they were added)
    reviewed = case_when(
      outcome == "Rare CNS (added after the release)" ~ "JLR",
      outcome != "different" ~ NA_character_,
      !is.na(reviewed_listed) ~ reviewed_listed,
      subtype_source %in% "clinical MYCN status" ~ "JLR",
      subtype_source %in% c("MYCN expression", "MYCN copy number") ~ "JLR",
      subtype_source %in% "manual review" ~ "JLR",
      TRUE ~ NA_character_),
    review_note = case_when(
      is.na(reviewed) ~ NA_character_,
      outcome == "Rare CNS (added after the release)" ~
        "Rare CNS subtype added after the release, subtyped and reviewed by the reviewer",
      !is.na(note_listed) ~ note_listed,
      subtype_source %in% "clinical MYCN status" ~ "MYCN subtype from the clinical report",
      subtype_source %in% c("MYCN expression", "MYCN copy number") ~
        "MYCN call changed by the NBL MYCN TPM fix (2026-10-06, 2026-10-08), after the release",
      subtype_source %in% "manual review" ~ "manual review (reviewer decision)",
      TRUE ~ NA_character_)
  ) %>%
  select(-reviewed_listed, -note_listed)

write_tsv(comparison, file.path(results_dir, "legacy_comparison.tsv"))

message("== Outcome overall")
print(count(comparison, outcome, sort = TRUE))

message("\n== Outcome by legacy group")
print(comparison %>%
        filter(!is.na(legacy_group)) %>%
        count(legacy_group, outcome) %>%
        pivot_wider(names_from = outcome, values_from = n, values_fill = 0),
      n = 40, width = Inf)

message("\n== Most common legacy -> new differences")
print(comparison %>% filter(outcome == "different") %>%
        count(legacy_subtype, new_subtype, sort = TRUE), n = 25)
