################################################################################
# mb-methylation-plots.R
#
# Plots of the medulloblastoma methylation arrays in a two-dimensional embedding
# (UMAP and PCA in 10-plot-mb-methylation.R),
# to check the subtypes of 05-validate-mb.R against the methylation profiles:
# all medulloblastoma, Group 3/4, WNT, and SHH with the marker genes (copy
# number and expression of MYCN, GLI2, CCND2, PTEN, TP53), TP53 status, age and
# the source of the SHH subtype (DKFZ v12.8, NIH v3.1, or age and molecular
# criteria).
#
# run_mb_embedding_plots(embed_fn, prefix, axes) reads the arrays (the parquet
# copy data/parquet/methyl-beta-values.parquet when it exists, otherwise the
# RDS in METHYL_BETA, default data/methyl-beta-values.rds) and writes
# plots/<prefix>_mb*.pdf and results/mb_<prefix>_coordinates.tsv.
#   embed_fn(mat, n_probes, seed): probes x arrays beta values -> matrix with two
#     columns (the coordinates) and the array IDs as row names
#   prefix: file name prefix, e.g. "umap"
#   axes: the two axis labels, e.g. c("UMAP1", "UMAP2"). embed_fn can set the
#     labels of a set of arrays (e.g. with the variance explained) in the
#     "axis_labels" attribute of the matrix it returns.
################################################################################

run_mb_embedding_plots <- function(embed_fn, prefix, axes) {
  axes_default <- axes
  suppressPackageStartupMessages({
    library(tidyverse)
    library(DBI)
  })

  root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
  module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
  results_dir <- file.path(module_dir, "results")
  input_dir <- file.path(module_dir, "input")
  plots_dir <- file.path(module_dir, "plots")
  dir.create(plots_dir, showWarnings = FALSE)
  source(file.path(module_dir, "utils", "duckdb-connect.R"))

  # the repo plot theme (white background)
  source(file.path(root_dir, "figures", "scripts", "figure_theme.R"))

  methyl_parquet <- file.path(root_dir, "data", "parquet", "methyl-beta-values.parquet")
  methyl_file <- Sys.getenv("METHYL_BETA",
                            file.path(root_dir, "data", "methyl-beta-values.rds"))
  stopifnot(file.exists(methyl_parquet) || file.exists(methyl_file))

  ### Arrays and annotations -----------------------------------------------------

  hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                   guess_max = 100000, show_col_types = FALSE)

  mb_validation <- read_tsv(file.path(results_dir, "mb_validation.tsv"),
                            show_col_types = FALSE)
  # tumors with a medulloblastoma subtype in the module table: mb_validation also
  # holds tumors called MB only by the expression classifier or by methylation
  # that another source overrides (a few Ewing, ETMR and sarcoma tumors)
  mb_tumors <- read_tsv(file.path(results_dir, "molecular_subtype.tsv"),
                        show_col_types = FALSE, guess_max = 100000) %>%
    filter(str_detect(replace_na(molecular_subtype, ""), "^MB, "),
           match_id %in% mb_validation$match_id) %>%
    pull(match_id)
  mb_arrays <- hist %>%
    filter(experimental_strategy == "Methylation",
           match_id %in% mb_tumors) %>%
    distinct(Kids_First_Biospecimen_ID, match_id, sample_id)

  # DKFZ v12.8 subclass and score of each array
  dkfz <- read_tsv(
    file.path(input_dir, "bti_methylation_classifications.tsv"),
    col_select = c(Kids_First_Biospecimen_ID = Bioassay_ID,
                   dkfz_subclass = `dkfz_v12.8_subclass`,
                   dkfz_score = `dkfz_v12.8_subclass_maxscore`),
    show_col_types = FALSE)

  # the module group (MB, SHH, MB, Group3, ...)
  subtype <- read_tsv(file.path(results_dir, "molecular_subtype.tsv"),
                      show_col_types = FALSE, guess_max = 100000) %>%
    distinct(match_id, molecular_subtype)

  # SHH subtype, its source, age, copy number, expression and TP53 per tumor
  shh_levels <- c("MB, SHH alpha", "MB, SHH beta", "MB, SHH gamma",
                  "MB, SHH delta", "MB, SHH unk")
  shh <- read_tsv(file.path(results_dir, "mb_shh_subtypes_w_molecular_data.tsv"),
                  show_col_types = FALSE, guess_max = 10000) %>%
    group_by(match_id) %>%
    summarise(across(everything(), ~ first(.x)), .groups = "drop") %>%
    transmute(
      match_id,
      shh_subtype = factor(if_else(molecular_subtype == "MB, SHH", "MB, SHH unk",
                                   molecular_subtype), levels = shh_levels),
      shh_source = case_when(
        !is.na(shh_methylation_source) ~ shh_methylation_source,
        !is.na(final_shh_subtype) ~ "age and molecular criteria",
        TRUE ~ "no subtype"),
      age_at_diagnosis_years,
      age_range = factor(case_when(
        age_at_diagnosis_years < 5 ~ "(0,5)",
        age_at_diagnosis_years < 10 ~ "[5,10)",
        age_at_diagnosis_years < 15 ~ "[10,15)",
        TRUE ~ ">=15"), levels = c("(0,5)", "[5,10)", "[10,15)", ">=15")),
      # a gene with no gain or amplification record is copy neutral
      across(c(consensus_CN_MYCN, consensus_CN_GLI2, consensus_CN_CCND2,
               consensus_CN_PTEN), ~ replace_na(.x, "neutral")),
      MYCN_TPM_zscore, GLI2_TPM_zscore, CCND2_TPM_zscore, PTEN_TPM_zscore,
      TP53_TPM_zscore, tp53_status)

  plot_df <- mb_arrays %>%
    left_join(dkfz, by = "Kids_First_Biospecimen_ID") %>%
    left_join(subtype, by = "match_id") %>%
    left_join(shh, by = "match_id")

  # colors to match the SHH subtypes in Garcia-Lopez 2020
  shh_colors <- c("MB, SHH alpha" = "aquamarine3", "MB, SHH beta" = "goldenrod2",
                  "MB, SHH gamma" = "royalblue1", "MB, SHH delta" = "plum4",
                  "MB, SHH unk" = "gray")

  ### Beta values of the arrays --------------------------------------------------

  if (file.exists(methyl_parquet)) {
    message("Reading the arrays from parquet")
    con <- connect_duckdb()
    src <- sprintf("read_parquet('%s')", normalizePath(methyl_parquet))
    available <- dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM %s", src))$column_name
    ids <- intersect(mb_arrays$Kids_First_Biospecimen_ID, available)
    methyl <- dbGetQuery(con, sprintf(
      "SELECT %s FROM %s",
      paste(sprintf('"%s"', c("Probe_ID", ids)), collapse = ", "), src))
    dbDisconnect(con, shutdown = TRUE)
  } else {
    message("Reading ", methyl_file, " (no parquet copy)")
    methyl <- readRDS(methyl_file)
    methyl <- methyl[, colnames(methyl) %in%
                       c("Probe_ID", mb_arrays$Kids_First_Biospecimen_ID)]
  }
  # the matrix stacks the arrays of different array types, so a probe can have a
  # second row that is NA for these arrays: keep the rows with values
  methyl <- methyl %>%
    filter(if_any(-Probe_ID, ~ !is.na(.x))) %>%
    distinct(Probe_ID, .keep_all = TRUE) %>%
    column_to_rownames("Probe_ID") %>%
    as.matrix()

  # embedding of a set of arrays: embed_fn(mat, n_probes, seed) takes the probes x
  # arrays beta values and returns a matrix with the two coordinates of each array
  # (row names are the array IDs)
  embed_of <- function(array_ids, n_probes = 20000, seed = 2024) {
    mat <- methyl[, colnames(methyl) %in% array_ids, drop = FALSE]
    emb <- embed_fn(mat, n_probes, seed)
    # the axis labels of the plots of this set of arrays (the plots are saved
    # after the embedding of their set and before the next one)
    labels <- attr(emb, "axis_labels")
    axes <<- if (is.null(labels)) axes_default else labels
    tibble(Kids_First_Biospecimen_ID = rownames(emb),
           Dim1 = emb[, 1], Dim2 = emb[, 2]) %>%
      left_join(plot_df, by = "Kids_First_Biospecimen_ID")
  }

  # small legends, so the plot area is larger
  small_legend <- theme(legend.title = element_text(size = 11, face = "bold"),
                        legend.text = element_text(size = 10),
                        legend.key.size = unit(0.4, "cm"),
                        legend.spacing.y = unit(0.15, "cm"),
                        legend.margin = margin(3, 3, 3, 3))

  save_plot <- function(p, name, width, height) {
    ggsave(file.path(plots_dir, paste0(name, ".pdf")),
           p + labs(x = axes[1], y = axes[2]) + small_legend, width = width, height = height)
  }

  # group labels of all medulloblastoma arrays, with SHH split by subtype
  group_label <- function(df) {
    df %>% mutate(group = if_else(!is.na(shh_subtype), as.character(shh_subtype),
                                  molecular_subtype))
  }

  ### All medulloblastoma ----------------------------------------------------------

  all_emb <- embed_of(mb_arrays$Kids_First_Biospecimen_ID) %>% group_label()
  save_plot(
    ggplot(all_emb, aes(Dim1, Dim2, color = group)) +
      geom_point(alpha = 0.7) +
      labs(color = "molecular subtype") +
      theme_Publication(),
    paste0(prefix, "_mb"), 6.5, 4.5)

  ### Group 3 / Group 4 ------------------------------------------------------------

  g34_ids <- plot_df %>%
    filter(molecular_subtype %in% c("MB, Group3", "MB, Group4")) %>%
    pull(Kids_First_Biospecimen_ID)
  g34_emb <- embed_of(g34_ids)
  save_plot(
    g34_emb %>%
      filter(str_detect(replace_na(dkfz_subclass, ""), "^MB_G34")) %>%
      ggplot(aes(Dim1, Dim2, color = dkfz_subclass, shape = molecular_subtype)) +
      geom_point(alpha = 0.7) +
      labs(color = "methylation subtype (DKFZ v12.8)", shape = "molecular subtype") +
      theme_Publication(),
    paste0(prefix, "_mb_group34"), 6.5, 4)

  ### WNT --------------------------------------------------------------------------

  wnt_ids <- plot_df %>% filter(molecular_subtype == "MB, WNT") %>%
    pull(Kids_First_Biospecimen_ID)
  wnt_emb <- embed_of(wnt_ids)
  save_plot(
    ggplot(wnt_emb, aes(Dim1, Dim2, color = dkfz_subclass)) +
      geom_point(alpha = 0.7, size = 3) +
      labs(color = "methylation subtype (DKFZ v12.8)") +
      theme_Publication(),
    paste0(prefix, "_mb_wnt"), 5.5, 3.5)

  ### SHH --------------------------------------------------------------------------

  shh_ids <- plot_df %>%
    filter(str_detect(replace_na(molecular_subtype, ""), "^MB, SHH")) %>%
    pull(Kids_First_Biospecimen_ID)
  shh_emb <- embed_of(shh_ids) %>%
    mutate(shh_subtype = fct_na_value_to_level(shh_subtype, "MB, SHH unk"))

  # subtype, age range and source of the subtype
  save_plot(
    ggplot(shh_emb, aes(Dim1, Dim2, color = shh_subtype, size = age_range,
                         shape = shh_source)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = shh_colors) +
      labs(color = "SHH subgroup", size = "age range (years)",
           shape = "source of the subtype") +
      guides(color = guide_legend(order = 1), shape = guide_legend(order = 2),
             size = guide_legend(order = 3)) +
      theme_Publication(),
    paste0(prefix, "_mb_shh"), 8, 5)

  # source of the subtype alone, to see where the non-methylation calls fall
  save_plot(
    ggplot(shh_emb, aes(Dim1, Dim2, color = shh_source)) +
      geom_point(alpha = 0.8, size = 3) +
      labs(color = "source of the SHH subtype") +
      theme_Publication(),
    paste0(prefix, "_mb_shh_source"), 6.5, 4)

  # DKFZ v12.8 subclass and its score
  save_plot(
    shh_emb %>%
      filter(str_detect(replace_na(dkfz_subclass, ""), "^MB_SHH")) %>%
      ggplot(aes(Dim1, Dim2, color = dkfz_subclass, alpha = dkfz_score)) +
      geom_point(size = 3) +
      scale_color_manual(values = c(MB_SHH_1 = "goldenrod2", MB_SHH_2 = "royalblue1",
                                    MB_SHH_3 = "aquamarine3", MB_SHH_4 = "plum4")) +
      labs(color = "methylation subtype (DKFZ v12.8)",
           alpha = "methylation subtype score") +
      theme_Publication(),
    paste0(prefix, "_mb_shh_methylation_subtype"), 6.5, 4)

  # marker genes: copy number status of MYCN, GLI2, CCND2 and PTEN
  save_plot(
    shh_emb %>%
      select(shh_subtype, Dim1, Dim2, MYCN = consensus_CN_MYCN,
             GLI2 = consensus_CN_GLI2, CCND2 = consensus_CN_CCND2,
             PTEN = consensus_CN_PTEN) %>%
      pivot_longer(c(MYCN, GLI2, CCND2, PTEN), names_to = "gene_name",
                   values_to = "CN_status") %>%
      ggplot(aes(Dim1, Dim2, color = shh_subtype, shape = CN_status)) +
      geom_point(alpha = 0.7, size = 3) +
      scale_color_manual(values = shh_colors) +
      facet_wrap(~gene_name, nrow = 1) +
      labs(color = "SHH subgroup", shape = "CN status") +
      theme_Publication(),
    paste0(prefix, "_mb_shh_cn_status"), 12, 4)

  # marker genes: expression (z-score of the TPM within the SHH tumors)
  expression_long <- shh_emb %>%
    select(Dim1, Dim2, MYCN = MYCN_TPM_zscore, GLI2 = GLI2_TPM_zscore,
           CCND2 = CCND2_TPM_zscore, PTEN = PTEN_TPM_zscore,
           TP53 = TP53_TPM_zscore) %>%
    pivot_longer(c(MYCN, GLI2, CCND2, PTEN, TP53), names_to = "gene_name",
                 values_to = "TPM_zscore") %>%
    # all missing when there are too few SHH tumors with RNA (the testing data): numeric,
    # so that the continuous scale still applies
    mutate(TPM_zscore = as.numeric(TPM_zscore))
  save_plot(
    ggplot(expression_long, aes(Dim1, Dim2, color = TPM_zscore)) +
      geom_point(size = 2.5) +
      scale_color_gradient2(low = "blue", mid = "gray90", high = "red",
                            midpoint = 0, na.value = "gray60") +
      facet_wrap(~gene_name, nrow = 1) +
      labs(color = "TPM z-score") +
      theme_Publication(),
    paste0(prefix, "_mb_shh_expression"), 14, 4)

  # TP53 status
  save_plot(
    ggplot(shh_emb, aes(Dim1, Dim2, color = shh_subtype, shape = tp53_status)) +
      geom_point(alpha = 0.7, size = 3) +
      scale_color_manual(values = shh_colors) +
      labs(color = "SHH subgroup", shape = "TP53 status") +
      theme_Publication(),
    paste0(prefix, "_mb_shh_tp53_status"), 6, 4)

  # age at diagnosis
  save_plot(
    ggplot(shh_emb, aes(Dim1, Dim2, color = age_at_diagnosis_years,
                         shape = shh_subtype)) +
      geom_point(alpha = 0.8, size = 3) +
      scale_color_viridis_c(option = "plasma", na.value = "gray60") +
      labs(color = "age at diagnosis (years)", shape = "SHH subgroup") +
      theme_Publication(),
    paste0(prefix, "_mb_shh_age"), 6.5, 4)

  write_tsv(
    bind_rows(mutate(all_emb, set = "all medulloblastoma"),
              mutate(g34_emb, set = "Group 3/4"),
              mutate(wnt_emb, set = "WNT"),
              mutate(shh_emb, set = "SHH")) %>%
      select(set, Kids_First_Biospecimen_ID, match_id, sample_id, Dim1, Dim2,
             molecular_subtype, shh_subtype, shh_source, dkfz_subclass,
             dkfz_score, age_at_diagnosis_years) %>%
      rename(!!axes_default[1] := Dim1, !!axes_default[2] := Dim2) %>%
      mutate(shh_subtype = as.character(shh_subtype)),
    file.path(results_dir, paste0("mb_", prefix, "_coordinates.tsv")))
  message("Wrote the plots to ", plots_dir)
}
