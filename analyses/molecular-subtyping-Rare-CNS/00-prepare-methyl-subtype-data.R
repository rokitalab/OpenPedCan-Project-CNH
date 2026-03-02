# Purpose:
# Read the DKFZ methylation annotation table, keep high-confidence Rare CNS
# methylation subclasses (score >= 0.8), and write a filtered table for
# downstream molecular subtyping.

library(tidyverse)
library(rprojroot)

# project root
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

methyl_file <- file.path(
  root_dir,
  "analyses",
  "molecular-subtyping-Rare-CNS",
  "subset-files",
  "mnp_v12.5_annotation_with_OPC_subtype_2025.tsv"
)

subtype_map_file <- file.path(
  root_dir,
  "analyses",
  "molecular-subtyping-Rare-CNS",
  "subset-files",
  "rare_cns_subtype_map.tsv"
)

independent_methyl_file <- file.path(
  root_dir,
  "data",
  "independent-specimens.methyl.primary-plus.tsv"
)

histology_file <- file.path(
  root_dir,
  "data",
  "histologies.tsv"
)

output_file <- file.path(
  root_dir,
  "analyses",
  "molecular-subtyping-Rare-CNS",
  "results",
  "methyl_rare_cns_subtyping.tsv"
)
# Read inputs

methyl_df <- readr::read_tsv(methyl_file, show_col_types = FALSE)
subtype_map <- readr::read_tsv(subtype_map_file, show_col_types = FALSE)
independent_methyl <- readr::read_tsv(independent_methyl_file, show_col_types = FALSE)
histologies <- readr::read_tsv(histology_file, show_col_types = FALSE)


required_map_cols <- c(
  "Abbreviation_internal_NIH",
  "dkfz_Abbreviation",
  "OPC_molecular_subtype"
)

missing_map_cols <- setdiff(required_map_cols, colnames(subtype_map))
if (length(missing_map_cols) > 0) {
  stop(
    "Missing required columns in rare_cns_subtype_map.tsv: ",
    paste(missing_map_cols, collapse = ", ")
  )
}

required_methyl_cols <- c(
  "Kids_First_Biospecimen_ID",
  "dkfz_v12_methylation_subclass",
  "dkfz_v12_methylation_subclass_score"
)

missing_methyl_cols <- setdiff(required_methyl_cols, colnames(methyl_df))
if (length(missing_methyl_cols) > 0) {
  stop(
    "Missing required columns in methyl annotation file: ",
    paste(missing_methyl_cols, collapse = ", "),
    "\nCheck the real column names and update the script."
  )
}

# Filter to Rare CNS methylation classes with high-confidence scores

rare_cns_methyl <- methyl_df %>%
  rename(
    dkfz_Abbreviation = dkfz_v12_methylation_subclass
  ) %>%
  inner_join(subtype_map, by = "dkfz_Abbreviation") %>%
  filter(dkfz_v12_methylation_subclass_score >= 0.8)


if ("Kids_First_Biospecimen_ID" %in% colnames(independent_methyl)) {
  rare_cns_methyl <- rare_cns_methyl %>%
    semi_join(independent_methyl, by = "Kids_First_Biospecimen_ID")
}


# Add histology 


if ("Kids_First_Biospecimen_ID" %in% colnames(histologies)) {
  rare_cns_methyl <- rare_cns_methyl %>%
    left_join(histologies, by = "Kids_First_Biospecimen_ID")
}


# write output


rare_cns_methyl <- rare_cns_methyl %>%
  arrange(OPC_molecular_subtype, desc(dkfz_v12_methylation_subclass_score))

readr::write_tsv(rare_cns_methyl, output_file)

message("Wrote filtered Rare CNS methylation table to: ", output_file)
message("Number of retained samples: ", nrow(rare_cns_methyl))