#' Methylation arrays of a likely mismatch
#'
#' 03-validate-subtypes.R decides the methylation mismatches
#' (results/methylation_mismatch_calls.tsv: found automatically from the lesions
#' and the pathology diagnosis, or decided by the reviewer in
#' input/possible_methylation_mixups.tsv). The arrays of a `likely` mismatch are
#' probably from another specimen, so no step uses them for the subtype; they are
#' to be removed from the release. `exclude_mismatched_arrays()` sets their
#' methyl_confidence to "excluded (likely mismatch)" in a table of calls.
likely_mismatch_arrays <- function(results_dir) {
  file <- file.path(results_dir, "methylation_mismatch_calls.tsv")
  if (!file.exists(file)) return(character())
  x <- readr::read_tsv(file, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  x <- x[x$mismatch_call %in% "likely", , drop = FALSE]
  unique(unlist(strsplit(stats::na.omit(x$methylation_biospecimen_ids), "; ", fixed = TRUE)))
}

exclude_mismatched_arrays <- function(calls, results_dir) {
  ids <- likely_mismatch_arrays(results_dir)
  calls$methyl_confidence[calls$Kids_First_Biospecimen_ID %in% ids] <- "excluded (likely mismatch)"
  calls
}
