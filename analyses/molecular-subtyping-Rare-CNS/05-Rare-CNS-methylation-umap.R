# Aylar Babaei for OpenPedCan 2026
#
# In this script, we generate a methylation UMAP for Rare CNS samples using
# high-confidence Rare CNS methylation subtype annotations.
#
# USAGE: Rscript --vanilla 05-Rare-CNS-methylation-umap.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(umap)
  library(ggplot2)
})

# Look for git root folder
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Paths
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "molecular-subtyping-Rare-CNS")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
subset_dir <- file.path(analysis_dir, "rare-cns-subset")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Files
hist_file <- file.path(data_dir, "histologies-base.tsv")
subtype_file <- file.path(results_dir, "rare_cns_subtyping.tsv")
methyl_file <- file.path(data_dir, "methyl-beta-values.rds")

# Read inputs
hist <- readr::read_tsv(hist_file, show_col_types = FALSE)
subtypes <- readr::read_tsv(subtype_file, show_col_types = FALSE)

# Keep methylation samples with Rare CNS subtype calls
hist_methyl <- hist %>%
  dplyr::filter(
    experimental_strategy == "Methylation",
    sample_type == "Tumor"
  ) %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    dkfz_v12_methylation_subclass,
    dkfz_v12_methylation_subclass_score,
    pathology_diagnosis
  ) %>%
  dplyr::distinct()

rare_cns_umap_samples <- subtypes %>%
  dplyr::filter(
    !is.na(molecular_subtype),
    molecular_subtype != "Rare CNS, To be classified"
  ) %>%
  dplyr::inner_join(
    hist_methyl,
    by = "Kids_First_Biospecimen_ID"
  ) %>%
  dplyr::distinct(Kids_First_Biospecimen_ID, .keep_all = TRUE)

message("Rare CNS methylation samples used for UMAP:")
print(rare_cns_umap_samples %>% dplyr::count(molecular_subtype, sort = TRUE))

# Load methylation beta matrix
methyl <- readRDS(methyl_file)

# Keep Probe_ID plus selected Rare CNS methylation biospecimens
rare_cns_methyl <- methyl[, colnames(methyl) %in% c("Probe_ID", rare_cns_umap_samples$Kids_First_Biospecimen_ID)]

rare_cns_methyl <- rare_cns_methyl %>%
  dplyr::distinct(Probe_ID, .keep_all = TRUE) %>%
  tibble::column_to_rownames("Probe_ID")

# Remove probes with zero variance or all NA
methyl_var <- apply(rare_cns_methyl, 1, var, na.rm = TRUE)
methyl_var[is.na(methyl_var)] <- 0

# Select most variable probes
n_var_probes <- min(20000, sum(methyl_var > 0))
var_probes <- names(sort(methyl_var, decreasing = TRUE)[seq_len(n_var_probes)])

message("Number of variable probes used: ", n_var_probes)

# Run UMAP
set.seed(2026)

n_samples <- ncol(rare_cns_methyl[var_probes, , drop = FALSE])
n_neighbors_use <- max(2, min(10, n_samples - 1))

umap_results <- umap::umap(
  t(rare_cns_methyl[var_probes, , drop = FALSE]),
  n_neighbors = n_neighbors_use
)

umap_plot_df <- data.frame(umap_results$layout) %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::left_join(
    rare_cns_umap_samples,
    by = "Kids_First_Biospecimen_ID"
  )

# Save coordinates
readr::write_tsv(
  umap_plot_df,
  file.path(results_dir, "rare_cns_umap_coordinates.tsv")
)

# Plot
rare_cns_umap_plot <- ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    fill = molecular_subtype
  )
) +
  geom_point(
    alpha = 0.85,
    size = 3.5,
    shape = 21,
    stroke = 0.8,
    color = "black"
  ) +
  labs(
    title = "Rare CNS methylation UMAP",
    x = "UMAP1",
    y = "UMAP2",
    fill = "Molecular subtype"
  ) +
  theme_bw()

ggsave(
  file.path(plots_dir, "rare_cns_methyl_umap.pdf"),
  rare_cns_umap_plot,
  width = 8,
  height = 5
)

ggsave(
  file.path(plots_dir, "rare_cns_methyl_umap.png"),
  rare_cns_umap_plot,
  width = 8,
  height = 5,
  dpi = 300
)

message("UMAP plot saved.")
message("Done.")