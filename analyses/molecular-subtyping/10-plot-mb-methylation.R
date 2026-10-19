################################################################################
# 10-plot-mb-methylation.R
#
# Check the medulloblastoma subtypes against the methylation profiles: UMAPs and
# PCAs of the most variable probes of the medulloblastoma arrays (all, Group 3/4,
# SHH, WNT), colored by the subtype of 05-validate-mb.R, and the SHH plots by the
# marker genes (copy number status and expression of MYCN, GLI2, CCND2, PTEN,
# TP53), TP53 status, age and the source of the SHH subtype (DKFZ v12.8, NIH
# v3.1, or age and molecular criteria). Tumors whose SHH subtype does not come
# from methylation should sit with the methylation-called tumors of the same
# subtype if the criteria are sound. The plots are in utils/mb-methylation-plots.R.
#
# UMAP: the 20,000 most variable probes (seeded, so the layout is the same at every
# run). PCA: the 20,000 most variable probes among the probes with a value in every
# array of the set, on the centered beta values; the axis labels give the variance
# explained by PC1 and PC2 of each set of arrays. PCA does not distort distances as
# UMAP can, so the subtypes should separate in it too.
#
# Plots (plots/): umap_mb* and pca_mb*: mb, mb_group34, mb_wnt, mb_shh,
# mb_shh_source, mb_shh_methylation_subtype, mb_shh_cn_status, mb_shh_expression,
# mb_shh_tp53_status, mb_shh_age.
# Coordinates: results/mb_umap_coordinates.tsv, results/mb_pca_coordinates.tsv.
#
# Needs the beta value matrix: only the medulloblastoma arrays are read from
# data/parquet/methyl-beta-values.parquet (bash scripts/download-parquet.sh
# methyl-beta-values), and otherwise from methyl-beta-values.rds (about 16 GB; set
# METHYL_BETA to its path, default data/methyl-beta-values.rds).
#
# Usage: Rscript --vanilla 10-plot-mb-methylation.R
################################################################################

suppressPackageStartupMessages(library(umap))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "analyses", "molecular-subtyping", "utils",
                 "mb-methylation-plots.R"))

# UMAP of the most variable probes (probes with an NA in an array are kept, as
# umap handles them)
embed_umap <- function(mat, n_probes, seed) {
  probes <- names(sort(apply(mat, 1, var, na.rm = TRUE), decreasing = TRUE))
  probes <- head(probes, n_probes)
  # seeded twice (R's generator and the umap configuration): the layout is
  # the same at every run
  set.seed(seed)
  umap::umap(t(mat[probes, ]), random_state = seed)$layout
}

embed_pca <- function(mat, n_probes, seed) {
  # prcomp needs every value: arrays of a different array type (for example 450k
  # next to EPIC) are missing some probes, so use the probes all arrays have
  complete <- mat[complete.cases(mat), , drop = FALSE]
  stopifnot(nrow(complete) >= 1000, ncol(complete) >= 3)
  variance <- apply(complete, 1, var)
  top <- complete[order(variance, decreasing = TRUE)[
    seq_len(min(n_probes, nrow(complete)))], , drop = FALSE]
  pca <- prcomp(t(top), center = TRUE, scale. = FALSE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  message(sprintf("PCA of %d arrays, %d probes (of %d complete): PC1 %.1f%%, PC2 %.1f%%",
                  ncol(top), nrow(top), nrow(complete), explained[1], explained[2]))
  scores <- pca$x[, 1:2, drop = FALSE]
  # axis labels with the variance explained
  attr(scores, "axis_labels") <- sprintf("PC%d (%.1f%%)", 1:2, explained[1:2])
  scores
}

run_mb_embedding_plots(embed_umap, prefix = "umap", axes = c("UMAP1", "UMAP2"))
run_mb_embedding_plots(embed_pca, prefix = "pca", axes = c("PC1", "PC2"))
