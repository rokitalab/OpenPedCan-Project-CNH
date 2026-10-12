################################################################################
# 07-compare-to-legacy.R
#
# Compare molecular_subtype_cns.tsv (06-compile-subtypes.R) with the legacy
# per-module subtypes compiled in
# molecular-subtyping-pathology/results/compiled_molecular_subtypes.tsv, per
# match_id. The legacy ", TP53" suffix is dropped, since TP53 is now an
# annotation, and DMG labels with `-mutant`, `-altered`, `EGFR` or `EZHIP` terms and `DHG, H3 G35-mutant` are compared as the legacy `DMG, H3 K28` and `DHG, H3 G35`.
# A legacy lesion-level label that equals `lesion_refinement` counts
# as agreement with the methylation family.
#
# Usage: Rscript --vanilla 07-compare-to-legacy.R
################################################################################

suppressPackageStartupMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
results_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS",
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
  file.path(root_dir, "analyses", "molecular-subtyping-pathology", "results",
            "compiled_molecular_subtypes.tsv"),
  col_types = cols(.default = col_character())
) %>%
  filter(!is.na(match_id), !is.na(molecular_subtype)) %>%
  transmute(match_id,
            legacy_subtype = str_remove(molecular_subtype, ", TP53$") %>%
              # labels renamed here to carry an HGG / LGG / GNT prefix
              recode(!!!legacy_renames),
            legacy_group = str_extract(legacy_subtype, "^[^,]+")) %>%
  distinct()

new <- read_tsv(file.path(results_dir, "molecular_subtype_cns.tsv"),
                guess_max = 10000, col_types = cols(.default = col_guess(),
                  tp53_snv_evidence = col_character(),
                  tp53_cnv_evidence = col_character(),
                  egfr_variants = col_character())) %>%
  # the legacy DMG, H3 K28 and DHG, H3 G35 cover both the methylation-only
  # label and the "-mutant" label (variant confirmed) used here
  mutate(new_subtype = molecular_subtype_cns %>%
           str_replace("^DMG, H3 K28(-mutant|-altered)?(, EGFR)?(, EZHIP)?$",
                       "DMG, H3 K28") %>%
           str_remove("(?<=H3 G35)-mutant$")) %>%
  select(match_id, new_subtype, subtype_source, lesion_refinement)

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
      legacy_subtype %in% lesion_refinement ~ "same (lesion refinement)",
      !legacy_assigned ~ "newly assigned",
      TRUE ~ "different"
    )
  )

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
