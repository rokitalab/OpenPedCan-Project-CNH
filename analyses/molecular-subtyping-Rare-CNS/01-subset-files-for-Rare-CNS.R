# Aylar Babaei for OpenPedCan 2026
#
# In this script, we identify candidate Rare CNS tumor biospecimens and gather
# high-confidence methylation subtype annotations for downstream Rare CNS
# subtyping.
#
# The Rare CNS module is primarily driven by high-confidence methylation
# subclassifications with dkfz_v12_methylation_subclass_score >= 0.8, with
# hallmark lesion evidence added in downstream scripts.
#
# USAGE: Rscript --vanilla 01-subset-files-for-Rare-CNS.R

library(tidyverse)

# Look for git root folder
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Get subset folder
subset_dir <- file.path(
  root_dir,
  "analyses",
  "molecular-subtyping-Rare-CNS",
  "rare-cns-subset"
)

# Create if it does not exist
if (!dir.exists(subset_dir)) {
  dir.create(subset_dir, recursive = TRUE)
}

# File from 00-Rare-CNS-select-pathology-dx.R that is used for pathology
# diagnosis inclusion/exclusion criteria and Rare CNS methylation mapping
path_dx_list <- jsonlite::fromJSON(
  file.path(subset_dir, "rare_cns_subtyping_path_dx_strings.json")
)

# Clinical file
clinical <- readr::read_tsv(
  file.path(root_dir, "data", "histologies-base.tsv"),
  guess_max = 100000
)

# Candidate Rare CNS tumor specimens based on pathology diagnosis
rare_cns_specimens_df <- clinical %>%
  dplyr::filter(
    stringr::str_detect(
      stringr::str_to_lower(pathology_diagnosis),
      paste(path_dx_list$include_path_dx, collapse = "|")
    ),
    stringr::str_detect(
      stringr::str_to_lower(pathology_diagnosis),
      paste(path_dx_list$exclude_path_dx, collapse = "|"),
      negate = TRUE
    ),
    sample_type == "Tumor"
  ) %>%
  dplyr::distinct()

# Write intermediate metadata file for inspection and downstream use
readr::write_tsv(
  rare_cns_specimens_df,
  file.path(subset_dir, "rare_cns_metadata.tsv")
)

# Assay-specific candidate biospecimens
rare_cns_dna_df <- rare_cns_specimens_df %>%
  dplyr::filter(
    experimental_strategy %in% c("WGS", "WXS", "Targeted Sequencing"),
    is.na(RNA_library)
  ) %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    Kids_First_Participant_ID,
    sample_id,
    match_id
  ) %>%
  dplyr::distinct()

rare_cns_rna_df <- rare_cns_specimens_df %>%
  dplyr::filter(
    experimental_strategy == "RNA-Seq" |
      (experimental_strategy == "Targeted Sequencing" & !is.na(RNA_library))
  ) %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    Kids_First_Participant_ID,
    sample_id,
    match_id
  ) %>%
  dplyr::distinct()

rare_cns_methyl_df <- clinical %>%
  dplyr::filter(
    experimental_strategy == "Methylation",
    sample_type == "Tumor"
  ) %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    Kids_First_Participant_ID,
    sample_id,
    match_id,
    dkfz_v12_methylation_subclass,
    dkfz_v12_methylation_subclass_score
  ) %>%
  dplyr::distinct()

readr::write_tsv(
  rare_cns_dna_df,
  file.path(subset_dir, "RareCNS_dna_biospecimen.tsv")
)

readr::write_tsv(
  rare_cns_rna_df,
  file.path(subset_dir, "RareCNS_rna_biospecimen.tsv")
)

readr::write_tsv(
  rare_cns_methyl_df,
  file.path(subset_dir, "RareCNS_methyl_biospecimen.tsv")
)

# High-confidence methylation subtypes
methyl_subtype_map <- path_dx_list$methyl_subtype_map %>%
  tibble::as_tibble()

rare_cns_methyl_subtyped <- clinical %>%
  dplyr::filter(
    experimental_strategy == "Methylation",
    sample_type == "Tumor",
    dkfz_v12_methylation_subclass_score >= 0.8
  ) %>%
  dplyr::left_join(
    methyl_subtype_map,
    by = c("dkfz_v12_methylation_subclass" = "dkfz_Abbreviation")
  ) %>%
  dplyr::filter(!is.na(OPC_molecular_subtype)) %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    Kids_First_Participant_ID,
    sample_id,
    match_id,
    pathology_diagnosis,
    pathology_free_text_diagnosis,
    dkfz_v12_methylation_subclass,
    dkfz_v12_methylation_subclass_score,
    Abbreviation_internal_NIH,
    OPC_molecular_subtype
  ) %>%
  dplyr::rename(
    Kids_First_Biospecimen_ID_Methyl = Kids_First_Biospecimen_ID,
    molecular_subtype_methyl = OPC_molecular_subtype
  ) %>%
  dplyr::distinct()

readr::write_tsv(
  rare_cns_methyl_subtyped,
  file.path(subset_dir, "RareCNS_methyl_subtypes.tsv")
)

# Optional console summaries
message("Candidate Rare CNS specimens by experimental strategy:")
print(rare_cns_specimens_df %>% dplyr::count(experimental_strategy))

message("Rare CNS methylation subtype counts:")
print(rare_cns_methyl_subtyped %>% dplyr::count(molecular_subtype_methyl, sort = TRUE))

message("Done.")