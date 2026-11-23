################################################################################
# 04-cell-line-fidelity.R
#
# Fidelity of the derived cell lines to their parental tumors, to track where a
# line departs from the tumor it came from. One row per cell line (match_id):
# the culture (composition, medium, passage), the line subtype and integrated
# diagnosis, the integrated diagnosis and subtype of the parental tumor (the
# solid-tissue tumors of the same sample, or failing that of the same
# participant), whether they agree, the confirming alterations of the tumor that
# the line lacks (from molecular-subtyping), the alterations seen in the line
# and the alterations that sibling lines of the same sample carry and this line
# does not (for instance a gene lost in one medium only).
#
# Usage: Rscript --vanilla 04-cell-line-fidelity.R
# Output: results/cell_line_fidelity.tsv
################################################################################

suppressMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
int_dir <- file.path(root_dir, "analyses", "molecular-subtyping-integrate", "results")
sub_dir <- file.path(root_dir, "analyses", "molecular-subtyping", "results")

collapse <- function(x) {
  x <- sort(unique(na.omit(x)))
  if (length(x) == 0) NA_character_ else paste(x, collapse = " | ")
}

hist <- read_tsv(file.path(int_dir, "histologies.tsv"), guess_max = 100000,
                 show_col_types = FALSE) %>%
  filter(sample_type == "Tumor", !is.na(match_id)) %>%
  select(match_id, Kids_First_Participant_ID, sample_id, cohort, composition,
         any_of(c("cell_line_composition", "cell_line_passage")), tumor_descriptor,
         molecular_subtype, integrated_diagnosis, cancer_group)

subtype <- read_tsv(file.path(sub_dir, "molecular_subtype.tsv"),
                    show_col_types = FALSE, col_types = cols(.default = "c")) %>%
  select(match_id, tumor_fidelity, tumor_alterations_shared,
         subtype_alterations_confirmed)
line_alt <- read_tsv(file.path(sub_dir, "cell_line_alterations.tsv"),
                     show_col_types = FALSE, col_types = cols(.default = "c")) %>%
  select(match_id, line_alterations)

is_line <- hist$composition %in% "Derived Cell Line"

# PBTA lines only: the neuroblastoma lines have no primary tumor to compare with
lines <- hist[is_line & hist$cohort %in% "PBTA", ] %>%
  group_by(match_id) %>%
  summarise(across(c(Kids_First_Participant_ID, sample_id, cohort, tumor_descriptor,
                     molecular_subtype, integrated_diagnosis, cancer_group),
                   first),
            across(any_of(c("cell_line_composition", "cell_line_passage")), collapse),
            .groups = "drop")

tumors <- hist[!is_line, ]
by_sample <- tumors %>% group_by(sample_id) %>%
  summarise(parental_subtype = collapse(molecular_subtype),
            parental_integrated_diagnosis = collapse(integrated_diagnosis),
            parental_cancer_group = collapse(cancer_group), .groups = "drop")
by_participant <- tumors %>% group_by(Kids_First_Participant_ID) %>%
  summarise(subtype_p = collapse(molecular_subtype),
            dx_p = collapse(integrated_diagnosis),
            group_p = collapse(cancer_group), .groups = "drop")

fidelity <- lines %>%
  left_join(by_sample, by = "sample_id", na_matches = "never") %>%
  left_join(by_participant, by = "Kids_First_Participant_ID", na_matches = "never") %>%
  mutate(parental_source = case_when(!is.na(parental_subtype) |
                                       !is.na(parental_integrated_diagnosis) ~ "same sample",
                                     !is.na(dx_p) | !is.na(subtype_p) ~ "same participant",
                                     TRUE ~ "no parental tumor"),
         parental_subtype = coalesce(parental_subtype, subtype_p),
         parental_integrated_diagnosis = coalesce(parental_integrated_diagnosis, dx_p),
         parental_cancer_group = coalesce(parental_cancer_group, group_p)) %>%
  select(-subtype_p, -dx_p, -group_p) %>%
  left_join(subtype, by = "match_id") %>%
  left_join(line_alt, by = "match_id") %>%
  mutate(
    integrated_diagnosis_agrees = case_when(
      is.na(parental_integrated_diagnosis) | is.na(integrated_diagnosis) ~ NA,
      TRUE ~ integrated_diagnosis == parental_integrated_diagnosis),
    tumor_alterations_not_in_line = if_else(
      str_detect(replace_na(tumor_alterations_shared, ""), "not found in line"),
      str_extract(tumor_alterations_shared, "not found in line: [^;]*") %>%
        str_remove("not found in line: "), NA_character_))

# alterations that sibling lines of the same sample carry and this line does not
alt_long <- fidelity %>%
  select(match_id, sample_id, line_alterations) %>%
  separate_rows(line_alterations, sep = "; ") %>%
  filter(!is.na(line_alterations), line_alterations != "") %>%
  mutate(alteration = line_alterations)
sibling <- fidelity %>%
  select(match_id, sample_id) %>%
  inner_join(alt_long %>% select(sample_id, alteration), by = "sample_id",
             relationship = "many-to-many") %>%
  distinct() %>%
  anti_join(alt_long %>% select(match_id, alteration), by = c("match_id", "alteration")) %>%
  group_by(match_id) %>%
  summarise(alterations_in_sibling_lines_only = paste(sort(unique(alteration)), collapse = "; "),
            .groups = "drop")

fidelity <- fidelity %>%
  left_join(sibling, by = "match_id") %>%
  mutate(fidelity = case_when(
    parental_source == "no parental tumor" ~ "no parental tumor",
    # the line has no subtype of its own beyond the default of its cancer group
    # (no methylation, or no lesion of the tumor found in the line)
    str_detect(replace_na(molecular_subtype, "To be classified"),
               "To be classified|H3 wildtype$|, NOS$") &
      !is.na(parental_subtype) &
      !str_detect(parental_subtype, "To be classified|H3 wildtype$|, NOS$") ~
      "line has a default subtype, the tumor a specific one",
    tumor_fidelity %in% "different" |
      integrated_diagnosis_agrees %in% FALSE ~ "differs from the parental tumor",
    !is.na(tumor_alterations_not_in_line) ~ "lacks a tumor alteration",
    !is.na(alterations_in_sibling_lines_only) ~ "differs from sibling lines",
    tumor_fidelity %in% c("same", "same group") ~ "same as the parental tumor",
    TRUE ~ "not assessable")) %>%
  arrange(fidelity, sample_id, cell_line_composition, cell_line_passage)

write_tsv(fidelity, file.path(int_dir, "cell_line_fidelity.tsv"))
message("Cell lines: ", nrow(fidelity))
print(count(fidelity, cohort, fidelity), n = Inf)
