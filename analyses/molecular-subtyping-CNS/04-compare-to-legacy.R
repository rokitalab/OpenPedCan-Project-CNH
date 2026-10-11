################################################################################
# 04-compare-to-legacy.R
#
# Compare subtypes from this module (methylation + lesion fallback) with the
# legacy per-module subtypes compiled in
# molecular-subtyping-pathology/results/compiled_molecular_subtypes.tsv, per
# match_id. The comparison is what justifies archiving the legacy modules.
#
# The legacy ", TP53" suffix is dropped before comparing, since TP53 is now an
# annotation rather than part of the subtype.
#
# Usage: Rscript --vanilla 04-compare-to-legacy.R
################################################################################

suppressPackageStartupMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS")
results_dir <- file.path(module_dir, "results")

# legacy labels that mean "no specific subtype was assigned"
unassigned_pattern <- "To be classified|, NOS$|wildtype$"

legacy <- read_tsv(
  file.path(root_dir, "analyses", "molecular-subtyping-pathology", "results",
            "compiled_molecular_subtypes.tsv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(match_id), !is.na(molecular_subtype)) %>%
  transmute(match_id,
            legacy_subtype = str_remove(molecular_subtype, ", TP53$"),
            legacy_group = str_extract(legacy_subtype, "^[^,]+")) %>%
  distinct()

# a match_id with more than one legacy subtype is reported, not silently resolved
legacy_multi <- legacy %>% count(match_id) %>% filter(n > 1)

validated <- read_tsv(file.path(results_dir, "methylation_subtype_validation.tsv"),
                      show_col_types = FALSE) %>%
  filter(!is.na(match_id)) %>%
  distinct(match_id, new_subtype = OPC_molecular_subtype) %>%
  mutate(new_source = "methylation")

fallback <- read_tsv(file.path(results_dir, "lesion_fallback_subtypes.tsv"),
                     show_col_types = FALSE) %>%
  filter(assign) %>%
  distinct(match_id, new_subtype = OPC_molecular_subtype) %>%
  mutate(new_source = "lesion fallback") %>%
  filter(!match_id %in% validated$match_id)

new <- bind_rows(validated, fallback) %>%
  group_by(match_id) %>%
  summarise(new_subtype = paste(sort(unique(new_subtype)), collapse = " | "),
            new_source = first(new_source), .groups = "drop")

comparison <- legacy %>%
  full_join(new, by = "match_id") %>%
  mutate(
    legacy_assigned = !is.na(legacy_subtype) &
      !str_detect(legacy_subtype, unassigned_pattern),
    outcome = case_when(
      is.na(legacy_subtype) ~ "new only",
      is.na(new_subtype) & legacy_assigned ~ "legacy only (specific)",
      is.na(new_subtype) ~ "legacy only (unassigned)",
      legacy_subtype == new_subtype ~ "same",
      !legacy_assigned ~ "newly assigned",
      TRUE ~ "different"
    )
  )

write_tsv(comparison, file.path(results_dir, "legacy_comparison.tsv"))

message("legacy match_ids with more than one subtype: ", nrow(legacy_multi))
message("\n== Outcome overall")
print(comparison %>% count(outcome, sort = TRUE))

message("\n== Outcome by legacy group")
print(comparison %>%
        filter(!is.na(legacy_group)) %>%
        count(legacy_group, outcome) %>%
        pivot_wider(names_from = outcome, values_from = n, values_fill = 0),
      n = 40, width = Inf)

message("\n== Most common legacy -> new differences")
print(comparison %>% filter(outcome == "different") %>%
        count(legacy_subtype, new_subtype, sort = TRUE), n = 40)

message("\n== Most common legacy-only specific subtypes")
print(comparison %>% filter(outcome == "legacy only (specific)") %>%
        count(legacy_subtype, sort = TRUE), n = 25)
