# Plot groups for the figures, from figures/palettes/plot-mapping.tsv.
#
# The table maps a cohort, a broad histology and a cancer group (as in the
# integrated histologies table) to a plot group with its colors and abbreviation,
# so every figure groups the tumors the same way. The original rows are for the
# PBTA cohort. Rows with a `status` were added (figures/scripts/extend-plot-mapping.R)
# for combinations that the original table did not have and are to be reviewed.
# The DGD, TCGA, PPTC, Maris and CPTAC cohorts are not in the table.
#
# Usage:
#   source(file.path(root_dir, "figures", "scripts", "plot_group.R"))
#   hist <- add_plot_group(hist)   # adds plot_group, plot_group_hex,
#                                  # broad_histology_display, broad_histology_hex
#                                  # and cancer_group_abbreviation

read_plot_mapping <- function(root_dir = rprojroot::find_root(rprojroot::has_dir(".git"))) {
  readr::read_tsv(file.path(root_dir, "figures", "palettes", "plot-mapping.tsv"),
                  col_types = readr::cols(.default = readr::col_character()),
                  na = character())
}

#' Add the plot group and its colors to a table with broad_histology and cancer_group
#'
#' A missing cancer group is the string "NA" in the table. A pair that is not in
#' the table gets no plot group and a warning that lists it.
add_plot_group <- function(df, mapping = read_plot_mapping()) {
  stopifnot(all(c("cohort", "broad_histology", "cancer_group") %in% names(df)))
  keep <- c("cohort", "broad_histology", "cancer_group", "broad_histology_display",
            "plot_group", "broad_histology_hex", "cancer_group_abbreviation",
            "plot_group_hex")
  key <- dplyr::mutate(df, .cancer_group_key = dplyr::coalesce(cancer_group, "NA"))
  mapping <- dplyr::select(mapping, dplyr::all_of(keep)) %>%
    dplyr::rename(.cancer_group_key = cancer_group)
  out <- dplyr::left_join(
    dplyr::select(key, -dplyr::any_of(setdiff(keep, c("cohort", "broad_histology", "cancer_group")))),
    mapping, by = c("cohort", "broad_histology", ".cancer_group_key"))
  missing <- out %>%
    dplyr::filter(is.na(plot_group), !is.na(broad_histology)) %>%
    dplyr::distinct(cohort, broad_histology, .cancer_group_key)
  if (nrow(missing) > 0) {
    warning(nrow(missing), " cohort, broad histology and cancer group combinations have no plot group in ",
            "figures/palettes/plot-mapping.tsv: ",
            paste(paste(missing$cohort, missing$broad_histology, missing$.cancer_group_key, sep = " / "),
                  collapse = "; "), call. = FALSE)
  }
  dplyr::select(out, -.cancer_group_key)
}
