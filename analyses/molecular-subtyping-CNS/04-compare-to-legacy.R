################################################################################
# 04-compare-to-legacy.R
#
# Compare molecular_subtype_cns.tsv (05-compile-subtypes.R) with the legacy
# per-module subtypes compiled in
# molecular-subtyping-pathology/results/compiled_molecular_subtypes.tsv, per
# match_id. The legacy ", TP53" suffix is dropped, since TP53 is now an
# annotation. A legacy lesion-level label that equals `lesion_refinement` counts
# as agreement with the methylation family.
#
# Usage: Rscript --vanilla 04-compare-to-legacy.R
################################################################################

suppressPackageStartupMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
results_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS",
                         "results")

# legacy labels that mean "no specific subtype was assigned"
unassigned <- "To be classified|, NOS$|wildtype$"

legacy <- read_tsv(
  file.path(root_dir, "analyses", "molecular-subtyping-pathology", "results",
            "compiled_molecular_subtypes.tsv"),
  col_types = cols(.default = col_character())
) %>%
  filter(!is.na(match_id), !is.na(molecular_subtype)) %>%
  transmute(match_id,
            legacy_subtype = str_remove(molecular_subtype, ", TP53$"),
            legacy_group = str_extract(legacy_subtype, "^[^,]+")) %>%
  distinct()

new <- read_tsv(file.path(results_dir, "molecular_subtype_cns.tsv"),
                guess_max = 10000, col_types = cols(.default = col_guess(),
                  tp53_snv_evidence = col_character(),
                  tp53_cnv_evidence = col_character(),
                  EGFR_variants = col_character())) %>%
  select(match_id, new_subtype = molecular_subtype_cns, subtype_source,
         lesion_refinement)

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
