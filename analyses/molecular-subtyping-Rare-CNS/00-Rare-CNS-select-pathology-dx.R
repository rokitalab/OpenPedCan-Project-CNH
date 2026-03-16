# Aylar Babaei for OpenPedCan 2026
#
# In this script we compile pathology diagnosis and methylation subtype
# terms/strings used as part of inclusion criteria for Rare CNS subtyping.
#
# The Rare CNS module will primarily use high-confidence methylation
# subclassifications (dkfz_v12_methylation_subclass_score >= 0.8) together
# with hallmark lesions in downstream scripts.
#
# USAGE: Rscript --vanilla 00-Rare-CNS-select-pathology-dx.R

# Detect the ".git" folder -- this will be in the project root directory.

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

# Inclusion criteria based on pathology diagnosis.
# These are intentionally broad candidate Rare CNS-related diagnoses so that
# downstream scripts can further refine using methylation subclass and
# hallmark lesion information.
include_path_dx <- stringr::str_to_lower(
  c(
    "Embryonal tumor",
    "Ewing sarcoma",
    "Neuroblastoma",
    "Sarcoma",
    "Glial-neuronal tumor NOS",
    "Low-grade glioma/astrocytoma",
    "Ganglioglioma",
    "Ependymoma"
  )
)

# Exclusion criteria for pathology diagnosis.
# These are entities that may overlap histologically with Rare CNS candidates
# but are expected to be handled in their own disease modules unless they are
# specifically reassigned downstream by Rare CNS methylation/lesion logic.
exclude_path_dx <- stringr::str_to_lower(
  c(
    "Medulloblastoma",
    "Atypical Teratoid/Rhabdoid Tumor",
    "Craniopharyngioma",
    "Chordoma"
  )
)

# Optional recode criteria on the basis of pathology_free_text_diagnosis.
# These terms can be used downstream to help flag candidate rare CNS entities.
recode_path_free_text <- stringr::str_to_lower(
  c(
    "mn1",
    "bcor",
    "foxr2",
    "cic",
    "dicer1",
    "plagl1",
    "patz1",
    "mycn",
    "myod1",
    "ews"
  )
)

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
  "NB_MYCN",
  "NB_TMM_NEG",
  "NB_TMM_POS",
  "NET_CXXC5",
  "NET_PATZ1",
  "NET_PLAGL1_FUS",
  "RMS_MYOD1"
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
  "MB_G3_Like",
  "HGNET_PLAG",
  "EWS",
  "NB_1",
  "NB_2s",
  "NB_2wo2s",
  "HGNET_CXXC5",
  "HGNET_PATZ",
  "EPN_ST_1",
  "RMS_MYOD1"
)

# OpenPedCan molecular subtype labels corresponding to the requested classes.
opc_molecular_subtype <- c(
  "Rare CNS, MN1::BEND2",
  "Rare CNS, ANTCON",
  "Rare CNS, ARMS",
  "Rare CNS, BCOR-fused",
  "Rare CNS, BCOR ITD",
  "Rare CNS, FOXR2 NBL",
  "Rare CNS, CIC SARC",
  "Rare CNS, DICER1 SARC",
  "Rare CNS, ERMS",
  "Rare CNS, BRD4::LEUTX",
  "Rare CNS, PLAG AMP",
  "Rare CNS, EWS",
  "Rare CNS, MYCN NBL",
  "Rare CNS, TMM NEG NBL",
  "Rare CNS, TMM POS NBL",
  "Rare CNS, MN1::CXXC5",
  "Rare CNS, PATZ1-fused",
  "Rare CNS, PLAGL1-fused",
  "Rare CNS, MYOD1 RMS"
)

# Create a lookup table to use in downstream scripts
methyl_subtype_map <- tibble::tibble(
  Abbreviation_internal_NIH = include_abbreviation_internal_nih,
  dkfz_Abbreviation = include_dkfz_abbreviation,
  OPC_molecular_subtype = opc_molecular_subtype
)

# Create a list with the strings and mapping tables used downstream
terms_list <- list(
  include_path_dx = include_path_dx,
  exclude_path_dx = exclude_path_dx,
  recode_path_free_text = recode_path_free_text,
  methyl_subtype_map = methyl_subtype_map
)

# Write to file
writeLines(
  jsonlite::prettify(jsonlite::toJSON(terms_list, pretty = TRUE, auto_unbox = TRUE)),
  output_file
)