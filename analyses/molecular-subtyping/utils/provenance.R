#' Provenance helpers shared by the subtyping scripts
#'
#' Notes use a vocabulary of their own (the legacy Notes strings are not
#' reused); `subtype_source` is the coded origin of the subtype.

#' Read the manual review table
#'
#' input/manual_subtype_review.tsv holds subtypes set by a pathology report or a
#' clinical review, keyed on sample_id. They have the highest precedence; the
#' automated call is kept in `automated_subtype`.
#'
#' @param module_dir The module directory.
#' @return sample_id, manual_subtype, manual_source and manual_note.
read_manual_review <- function(module_dir) {
  manual_review <- read_tsv(file.path(module_dir, "input", "manual_subtype_review.tsv"),
                            show_col_types = FALSE) %>%
    transmute(sample_id, manual_subtype = molecular_subtype,
              manual_source = source, manual_note = note)
  stopifnot(!anyDuplicated(manual_review$sample_id))
  manual_review
}

#' match_id for each manually reviewed sample (tumor rows of histologies-base)
reviewed_tumor_ids <- function(manual_review, root_dir) {
  read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
           guess_max = 100000, show_col_types = FALSE) %>%
    filter(sample_type == "Tumor", !is.na(match_id),
           sample_id %in% manual_review$sample_id) %>%
    distinct(sample_id, match_id, Kids_First_Participant_ID) %>%
    group_by(sample_id) %>%
    slice(1) %>%
    ungroup()
}

#' Notes text from the source of a subtype
#'
#' @param subtype_source The coded source.
#' @param molecular_subtype The subtype label (used to tell ETMR from the other
#'   embryonal rules).
subtype_notes <- function(subtype_source, molecular_subtype) {
  case_when(
    subtype_source == "manual review" ~
      "Subtype from manual review of pathology report",
    subtype_source == "methylation" ~ "Subtype from methylation",
    subtype_source == "lesion overrides methylation" ~
      "Subtype from lesion, overriding methylation",
    subtype_source %in% c("lesion fallback", "LGG lesion", "SMARCB1 focal loss") ~
      "Subtype from lesion",
    subtype_source == "lesion negative" ~
      "Subtype by exclusion, no H3 or IDH variant",
    subtype_source == "expression classifier" ~
      "Subtype from expression classifier",
    subtype_source == "pathology free text" ~
      "Subtype from pathology free text diagnosis",
    subtype_source == "pathology diagnosis" ~
      "Subtype from pathology diagnosis and primary site",
    subtype_source == "EWSR1 fusion" ~ "Subtype from EWSR1 fusion",
    subtype_source == "MYCN copy number" ~ "Subtype from MYCN copy number",
    subtype_source == "MYCN expression" ~ "Subtype from MYCN expression",
    subtype_source == "clinical MYCN status" ~ "Subtype from clinical MYCN status",
    subtype_source == "embryonal rules" &
      str_detect(replace_na(molecular_subtype, ""), "^ETMR") ~
      "Subtype from fusion, copy number and expression",
    subtype_source == "embryonal rules" ~
      "Subtype by exclusion, embryonal diagnosis",
    TRUE ~ NA_character_
  )
}
