# Aylar Babaei, Jo Lynne Rokita for OpenPedCan 2026
#
# USAGE: Rscript --vanilla 00-Rare-CNS-select-tumors.R


root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# JSON file containing the terms/strings
output_file <- file.path(
  root_dir,
  "analyses",
  "molecular-subtyping-Rare-CNS",
  "rare-cns-subset",
  "rare_cns_subtyping_path_dx_strings.json"
)

# Create output directory if it does not exist
if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}

# Inclusion criteria for high-confidence DKFZ methylation subclasses
# requested in the issue.
include_dkfz_abbreviation <- c(
  "ABM_MN1",
  "ANTCON",
  "ARMS",
  "CNS_BCOR_FUS",
  "CNS_BCOR_ITD",
  "CNS_NB_FOXR2",
  "CNS_SARC_CIC",
  "CNS_SARC_DICER",
  "ERMS",
  "ET_BRD4_LEUTX",
  "ET_PLAG",
  "EWS",
  "NET_CXXC5",
  "NET_PATZ1",
  "NET_PLAGL1_FUS",
  "RMS_MYOD1",
  "CRINET"
)

# Companion NIH/internal abbreviations from the issue.
include_abbreviation_internal_nih <- c(
  "HGNET_MN1",
  "ANTCON",
  "RMS_alevolar",
  "EP300_BCOR",
  "HGNET_BCOR",
  "CNS_NB_FOXR2",
  "EFT_CIC",
  "CNS_SCS_DICER",
  "RMS_embryonal",
  "ET_BRD4_LEUTX",
  "HGNET_PLAG",
  "EWS",
  "HGNET_CXXC5",
  "HGNET_PATZ",
  "PLAGL1_FUS",
  "RMS_MYOD1",
  "CRINET"
)

exclude_diagnoses <- c("Meningioma", 
                "Choroid plexus carcinoma", 
                "Brainstem glioma- Diffuse intrinsic pontine glioma",
                "Medulloblastoma",
                "Chordoma")

exclude_subtypes <- c("DMG_K27")

# Create a list with the strings and mapping tables used downstream
terms_list <- list(
  include_dkfz_subtypes = include_dkfz_abbreviation,
  include_nih_subtypes = include_abbreviation_internal_nih,
  exclude_dx = exclude_diagnoses,
  exclude_sub = exclude_subtypes
  
)

# Write to file
writeLines(
  jsonlite::prettify(jsonlite::toJSON(terms_list, pretty = TRUE, auto_unbox = TRUE)),
  output_file
)