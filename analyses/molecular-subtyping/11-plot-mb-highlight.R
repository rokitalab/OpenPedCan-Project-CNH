################################################################################
# 11-plot-mb-highlight.R
#
# Highlight chosen medulloblastoma tumors in the UMAP and PCA of the methylation
# arrays (coordinates of
# 10-plot-mb-methylation.R), to see whether tumors whose subtype changed, or
# that the expression classifier calls another group, sit in the cluster of their
# methylation subtype. The default is the tumors whose Group 3/4 or SHH call
# differs from the previous release and from the expression classifier.
#
# Usage: Rscript --vanilla 11-plot-mb-highlight.R [sample_id ...]
# Output: plots/umap_mb_highlight.pdf, plots/pca_mb_highlight.pdf (all
#   medulloblastoma) and plots/umap_mb_group34_highlight.pdf,
#   plots/pca_mb_group34_highlight.pdf (Group 3/4)
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
source(file.path(root_dir, "figures", "scripts", "figure_theme.R"))

args <- commandArgs(trailingOnly = TRUE)
highlight <- if (length(args) > 0) args else
  c("7316-4956", "7316-358", "7316-6794", "7316-6806", "7316-1772")

validation <- read_tsv(file.path(module_dir, "results", "mb_validation.tsv"),
                       show_col_types = FALSE) %>%
  filter(sample_id %in% highlight) %>%
  transmute(sample_id,
            label = paste0(sample_id, " (", str_remove(mb_methylation_group, "MB, "),
                           " by methylation, ", str_remove(mb_expression_group, "MB, "),
                           " by expression)"))

group_colors <- c("MB, Group3" = "#E69F00", "MB, Group4" = "#56B4E9",
                  "MB, SHH" = "#D55E00", "MB, WNT" = "#009E73",
                  "MB, Medullomyoblastoma" = "grey40")

plot_highlight <- function(file, x, y, set_name, out) {
  d <- read_tsv(file.path(module_dir, "results", file), show_col_types = FALSE) %>%
    filter(set == set_name) %>%
    left_join(validation, by = "sample_id") %>%
    mutate(highlighted = !is.na(label))
  p <- ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_point(data = filter(d, !highlighted),
               aes(color = molecular_subtype), size = 1.3, alpha = 0.7) +
    geom_point(data = filter(d, highlighted), aes(fill = molecular_subtype),
               shape = 21, color = "black", size = 3.8, stroke = 1.1) +
    geom_text_repel(data = filter(d, highlighted), aes(label = label),
                    size = 2.6, min.segment.length = 0, box.padding = 0.6,
                    max.overlaps = Inf) +
    scale_color_manual(values = group_colors, name = "Subtype") +
    scale_fill_manual(values = group_colors, guide = "none") +
    theme_Publication(base_size = 11) +
    theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9))
  ggsave(file.path(module_dir, "plots", out), p, width = 8, height = 6.5)
  invisible(d)
}

plot_highlight("mb_umap_coordinates.tsv", "UMAP1", "UMAP2", "all medulloblastoma",
               "umap_mb_highlight.pdf")
plot_highlight("mb_umap_coordinates.tsv", "UMAP1", "UMAP2", "Group 3/4",
               "umap_mb_group34_highlight.pdf")
plot_highlight("mb_pca_coordinates.tsv", "PC1", "PC2", "all medulloblastoma",
               "pca_mb_highlight.pdf")
plot_highlight("mb_pca_coordinates.tsv", "PC1", "PC2", "Group 3/4",
               "pca_mb_group34_highlight.pdf")
message("Highlighted: ", paste(validation$sample_id, collapse = ", "))
