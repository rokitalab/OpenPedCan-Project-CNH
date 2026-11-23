################################################################################
# 02-compare-to-release.R
#
# List the molecular subtypes that changed between a previous release table and
# the current integrated table (results/histologies.tsv), one row per tumor
# (match_id), with the evidence behind the new subtype and a verdict on whether
# the change makes sense.
#
# Change type (the previous label is `old`, the current label `new`):
#   vocabulary         the same subtype under a new name (TP53 suffix removed,
#                      H3 K28 -> H3 K28-mutant, ...; input/subtype_vocabulary.tsv)
#   newly subtyped     old is missing or a default ("To be classified", "..., NOS",
#                      "..., wildtype"), new is a specific subtype
#   refined            old is a subtype, new is a more specific or different
#                      subtype of the same cancer group
#   reclassified       the cancer group changed
#   lost               new is missing or a default
#
# Verdict (rules, from the evidence columns of molecular-subtyping):
#   no change in meaning   vocabulary
#   supported              methylation with a concordant or compatible diagnosis,
#                          lesion confirmation, a manual review, or an MYCN call
#                          from copy number or the clinical report
#   supported, diagnosis   methylation or lesion call and the diagnosis cannot
#   uninformative          judge it (no informative diagnosis)
#   review                 the diagnosis is discordant, a lesion contradicts it,
#                          the classifiers disagree, or the call carries a
#                          review flag
#
# Usage: Rscript --vanilla 02-compare-to-release.R [release histologies file]
#   (default data/v15/histologies.tsv)
# Outputs (results/): subtype_changes_from_release.tsv (one row per tumor),
#   subtype_change_summary.tsv (change type x verdict),
#   subtype_change_transitions.tsv (old -> new label counts with the verdicts),
#   subtype_changes_needing_review.tsv (the tumors with the verdict "review"),
#   cancer_group_changes.tsv (biospecimens whose cancer group changed, by the
#   previous and the current cancer group and integrated diagnosis)
################################################################################

suppressMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-integrate")
results_dir <- file.path(module_dir, "results")

args <- commandArgs(trailingOnly = TRUE)
release_file <- if (length(args) >= 1) args[1] else
  file.path(root_dir, "data", "v15", "histologies.tsv")

# the label means "no specific subtype" (as in 09-compare-to-legacy.R, except
# that "H3 wildtype" is a subtype: assigned by exclusion of H3 and IDH variants)
unassigned_regex <- "To be classified|(?<!ETMR), NOS$|(?<!H3 )wildtype$"
is_unassigned <- function(x) is.na(x) | x == "NA" | str_detect(x, unassigned_regex)

# an empty replacement removes the match, and spaces in a pattern are kept
vocabulary <- read_tsv(file.path(module_dir, "input", "subtype_vocabulary.tsv"),
                       show_col_types = FALSE, col_types = cols(.default = "c"),
                       na = character(), trim_ws = FALSE)
normalize_label <- function(x) {
  for (i in seq_len(nrow(vocabulary))) {
    x <- str_replace(x, vocabulary$pattern[i], vocabulary$replacement[i])
  }
  x
}

### Release and current tables -----------------------------------------------------

keep <- c("Kids_First_Biospecimen_ID", "match_id", "sample_id", "sample_type",
          "cohort", "pathology_diagnosis", "molecular_subtype", "cancer_group",
          "integrated_diagnosis", "broad_histology", "short_histology")
release <- read_tsv(release_file, col_select = any_of(keep), guess_max = 100000,
                    show_col_types = FALSE)
current <- read_tsv(file.path(results_dir, "histologies.tsv"),
                    col_select = any_of(keep), guess_max = 100000,
                    show_col_types = FALSE)

collapse_labels <- function(x) {
  x <- sort(unique(na.omit(x)))
  if (length(x) == 0) NA_character_ else paste(x, collapse = " | ")
}

old <- release %>%
  filter(sample_type == "Tumor", !is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(old_subtype = collapse_labels(molecular_subtype),
            old_subtype_normalized = collapse_labels(normalize_label(molecular_subtype)),
            old_cancer_group = collapse_labels(cancer_group),
            old_integrated_diagnosis = collapse_labels(integrated_diagnosis),
            .groups = "drop")
new <- current %>%
  filter(sample_type == "Tumor", !is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(sample_id = first(sample_id), cohort = first(cohort),
            pathology_diagnosis = first(pathology_diagnosis),
            new_subtype = collapse_labels(molecular_subtype),
            new_subtype_normalized = collapse_labels(normalize_label(molecular_subtype)),
            new_cancer_group = collapse_labels(cancer_group),
            new_integrated_diagnosis = collapse_labels(integrated_diagnosis),
            new_broad_histology = collapse_labels(broad_histology),
            new_short_histology = collapse_labels(short_histology),
            .groups = "drop")

# evidence behind the new subtype
evidence_cols <- c("match_id", "subtype_source", "methyl_class", "methyl_score",
                   "methyl_conflict", "dx_concordance", "lesion_status",
                   "lesion_evidence", "lesion_contradiction", "subtype_review_flag",
                   "manual_conflict", "manual_source", "needs_review", "diagnosis_review",
                   "subtype_alterations_confirmed")
evidence <- read_tsv(
  file.path(root_dir, "analyses", "molecular-subtyping", "results",
            "molecular_subtype.tsv"),
  col_select = any_of(evidence_cols), guess_max = 100000, show_col_types = FALSE,
  col_types = cols(.default = col_guess(), methyl_conflict = col_character(),
                   needs_review = col_character(), manual_conflict = col_character(),
                   lesion_contradiction = col_character())) %>%
  distinct(match_id, .keep_all = TRUE)

changes <- inner_join(old, new, by = "match_id") %>%
  filter(replace_na(old_subtype, "") != replace_na(new_subtype, "")) %>%
  left_join(evidence, by = "match_id")

# one string per row from several columns, dropping the missing and the repeats
join_reasons <- function(...) {
  parts <- cbind(...)
  vapply(seq_len(nrow(parts)),
         function(i) {
           pieces <- str_trim(unlist(str_split(na.omit(parts[i, ]), ";")))
           paste(unique(pieces[pieces != ""]), collapse = "; ")
         }, "")
}

### Change type and verdict ----------------------------------------------------------

changes <- changes %>%
  mutate(
    old_unassigned = is_unassigned(old_subtype),
    new_unassigned = is_unassigned(new_subtype),
    group_changed = replace_na(old_cancer_group, "") != replace_na(new_cancer_group, ""),
    change_type = case_when(
      !is.na(old_subtype) & !is.na(new_subtype) &
        old_subtype_normalized == new_subtype_normalized ~ "vocabulary",
      old_unassigned & !new_unassigned ~ "newly subtyped",
      !old_unassigned & new_unassigned ~ "lost",
      group_changed ~ "reclassified",
      TRUE ~ "refined"
    ),
    strong_source = subtype_source %in% c("manual review", "clinical MYCN status",
                                          "MYCN copy number"),
    # a manual review is a decision: its difference from the automated call is
    # recorded (manual_note) and is not a reason to review it again
    flagged = !subtype_source %in% "manual review" &
      (!is.na(subtype_review_flag) | !is.na(lesion_contradiction) |
         (dx_concordance %in% "discordant" & is.na(diagnosis_review)) |
         methyl_conflict %in% "TRUE"),
    supported = strong_source | lesion_status %in% "confirmed" |
      !is.na(diagnosis_review) |
      (subtype_source %in% "methylation" &
         dx_concordance %in% c("concordant", "compatible")),
    verdict = case_when(
      change_type == "vocabulary" ~ "no change in meaning",
      flagged ~ "review",
      supported ~ "supported",
      TRUE ~ "supported, diagnosis uninformative"
    ),
    reason = case_when(
      change_type == "vocabulary" ~ "same subtype, new name",
      verdict == "supported" & !is.na(diagnosis_review) ~
        paste0("methylation ", methyl_class, ", discordant diagnosis reviewed: ",
               str_remove(diagnosis_review, ":.*")),
      verdict == "supported" & subtype_source %in% "methylation" &
        lesion_status %in% "confirmed" ~
        paste0("methylation ", methyl_class, " with lesion confirmation (",
               lesion_evidence, ")"),
      verdict == "supported" & subtype_source %in% "methylation" ~
        paste0("methylation ", methyl_class, " (", signif(methyl_score, 2),
               "), diagnosis ", dx_concordance),
      verdict == "supported" ~ paste0(subtype_source,
                                      if_else(lesion_status %in% "confirmed",
                                              paste0(": ", lesion_evidence), "")),
      verdict == "review" ~ join_reasons(
        if_else(dx_concordance %in% "discordant", "diagnosis discordant", NA_character_),
        subtype_review_flag, lesion_contradiction,
        if_else(methyl_conflict %in% "TRUE", "classifier versions disagree", NA_character_)),
      TRUE ~ paste0(subtype_source, ", diagnosis not informative")
    )
  ) %>%
  select(match_id, sample_id, cohort, pathology_diagnosis, old_subtype, new_subtype,
         change_type, verdict, reason, old_cancer_group, new_cancer_group,
         old_integrated_diagnosis, new_integrated_diagnosis, new_broad_histology,
         new_short_histology, subtype_source, methyl_class, methyl_score,
         dx_concordance, lesion_status, lesion_evidence, subtype_review_flag) %>%
  arrange(change_type, verdict, old_subtype, new_subtype)

write_tsv(changes, file.path(results_dir, "subtype_changes_from_release.tsv"))

write_tsv(changes %>% filter(verdict == "review"),
          file.path(results_dir, "subtype_changes_needing_review.tsv"))

# biospecimens whose cancer group changed
inner_join(
  release %>% select(Kids_First_Biospecimen_ID,
                     previous_cancer_group = cancer_group,
                     previous_integrated_diagnosis = integrated_diagnosis),
  current %>% select(Kids_First_Biospecimen_ID, cancer_group, integrated_diagnosis),
  by = "Kids_First_Biospecimen_ID"
) %>%
  filter(replace_na(previous_cancer_group, "") != replace_na(cancer_group, "")) %>%
  count(previous_cancer_group, cancer_group, previous_integrated_diagnosis,
        integrated_diagnosis, name = "biospecimens") %>%
  arrange(desc(biospecimens)) %>%
  write_tsv(file.path(results_dir, "cancer_group_changes.tsv"))

summary_table <- changes %>% count(change_type, verdict) %>%
  pivot_wider(names_from = verdict, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric))))
write_tsv(summary_table, file.path(results_dir, "subtype_change_summary.tsv"))

transitions <- changes %>%
  count(change_type, old_subtype, new_subtype, verdict) %>%
  pivot_wider(names_from = verdict, values_from = n, values_fill = 0) %>%
  mutate(tumors = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(tumors))
write_tsv(transitions, file.path(results_dir, "subtype_change_transitions.tsv"))

message("Release: ", release_file)
message("Tumors in both: ", nrow(inner_join(old, new, by = "match_id")),
        "; subtype changed: ", nrow(changes))
print(summary_table, width = Inf)
