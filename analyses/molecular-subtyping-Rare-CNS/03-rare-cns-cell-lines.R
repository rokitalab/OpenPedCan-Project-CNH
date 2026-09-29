#!/usr/bin/env Rscript

# Identify derived cell lines from participants/samples classified by the
# Rare-CNS molecular-subtyping workflow.

library(tidyverse)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-Rare-CNS")

rare_subtypes_file <- file.path(
  module_dir, "results", "rare-cns-molecular-subtypes.tsv"
)
histologies_file <- file.path(root_dir, "data", "histologies.tsv")
output_file <- file.path(
  module_dir, "results", "rare-cns-cell-line-histologies.tsv"
)

if (!file.exists(rare_subtypes_file)) {
  stop("Rare-CNS subtype results are missing: ", rare_subtypes_file)
}

if (!file.exists(histologies_file)) {
  stop("Final histologies file is missing: ", histologies_file)
}

# A participant can have multiple tumors/events. Match on sample_id as well as
# participant ID so a cell line is linked only to the Rare-CNS-classified tumor.
rare_subtypes <- readr::read_tsv(rare_subtypes_file, show_col_types = FALSE) %>%
  select(
    Kids_First_Participant_ID,
    sample_id,
    rare_cns_molecular_subtype = molecular_subtype
  ) %>%
  distinct()

cell_lines <- readr::read_tsv(histologies_file, show_col_types = FALSE) %>%
  filter(composition == "Derived Cell Line") %>%
  select(
    Kids_First_Participant_ID,
    Kids_First_Biospecimen_ID,
    sample_id,
    aliquot_id,
    experimental_strategy,
    composition,
    cell_line_composition,
    cell_line_passage,
    tumor_descriptor,
    pathology_diagnosis,
    histologies_molecular_subtype = molecular_subtype,
    cancer_group,
    match_id
  )

rare_cell_lines <- rare_subtypes %>%
  inner_join(cell_lines, by = c("Kids_First_Participant_ID", "sample_id")) %>%
  select(
    Kids_First_Participant_ID,
    Kids_First_Biospecimen_ID,
    sample_id,
    aliquot_id,
    experimental_strategy,
    composition,
    cell_line_composition,
    cell_line_passage,
    tumor_descriptor,
    pathology_diagnosis,
    histologies_molecular_subtype,
    cancer_group,
    match_id,
    rare_cns_molecular_subtype
  ) %>%
  arrange(Kids_First_Participant_ID, sample_id, match_id, experimental_strategy)

readr::write_tsv(rare_cell_lines, output_file, na = "NA")

message(
  "Wrote ", nrow(rare_cell_lines), " derived-cell-line histology rows to ",
  output_file
)
