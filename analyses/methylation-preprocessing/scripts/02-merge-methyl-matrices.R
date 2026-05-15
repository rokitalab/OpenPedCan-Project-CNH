# Merges methylation beta-values, m-values, and cp-values matrices for 
# all pre-processed array datasets.

# Jessica Daggett
# 05/12/2026

#you are very close with this - just need to update file naming conventions and add lines to 
#remove memory heavy objects that you do not need as you go. 



# set up optparse options
option_list <- list(
  make_option(opt_str = "--output_dir", type = "character", default = NULL,
              help = "The directory containing the output .parquet files",
              metavar = "character"),
  make_option(opt_str = "--output_prefix", type = "character", default = NULL,
              help = "output file prefix",
              metavar = "character"),
)

# establish base dir
#root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

out_dir <- 'test-out'
out_pref <- 'test'


#define the data types - these will be looped through later
data_types <- c('beta-values-masked', 'm-values-masked', 'm-values-unmasked', 'cn-values')


#find out what array types are present - gets files
files <- list.files(out_dir, full.names = FALSE)
array_types <- character()

if (any(grepl("EPICv2", files, ignore.case = TRUE))) {
  array_types <- c(array_types, "EPICv2")
}
if (
  any(grepl("EPICv1", files, ignore.case = TRUE)) ||
  any(grepl("EPIC", files, ignore.case = TRUE) &
      !grepl("EPICv2", files, ignore.case = TRUE))
) {
  array_types <- c(array_types, "EPICv1")
}
if (any(grepl("450k", files, ignore.case = TRUE))) {
  array_types <- c(array_types, "450k")
}

if ("EPICv2" %in% array_types) {
  # --- Load detP parquet ---
  detP_epicv2_fn <- file.path(
    out_dir,
    paste0(out_pref, "-IlluminaHumanMethylationEPICv2-methyl-p-values.parquet")
  )
  detP <- open_dataset(detP_epicv2_fn)
  # --- Identify assay columns ---
  cols <- setdiff(colnames(detP), "ProbeID")
  # --- Collect only required data ---
  detP_df <- detP %>%
    select(ProbeID, all_of(cols)) %>%
    collect()
  # --- Compute mean detection p-values (row-wise) ---
  detP_df$detP_mean <- rowMeans(detP_df[, cols], na.rm = TRUE)
  # Optional but recommended: drop wide columns to save memory
  detP_df[, cols] <- NULL
  
  # --- Select best probe per Probe_base ---
  epicv2unique <- detP_df %>%
    mutate(
      Probe_base = sub("_.*$", "", ProbeID)
    ) %>%
    select(ProbeID, Probe_base, detP_mean) %>%
    group_by(Probe_base) %>%
    slice_min(order_by = detP_mean, n = 1, with_ties = FALSE) %>%
    ungroup()
  rm(detP, detP_df)
  gc()
}

#start with one dataset - try 
if (all(c("EPICv1", "EPICv2") %in% array_types)) {
  
  # --- File paths ---
  epicv1_vals_fn <- file.path(
    out_dir,
    paste0(out_pref, "-IlluminaHumanMethylationEPIC-methyl-p-values.parquet")
  )
  # --- Open datasets (lazy) ---
  epicv1_vals <- open_dataset(epicv1_vals_fn)
  # --- Get EPICv1 probes (reference set) ---
  v1_probes <- epicv1_vals %>%
    select(ProbeID) %>%
    distinct() %>%
    collect() %>%
    pull(ProbeID)
  # --- Find common probes (via Probe_base mapping) ---
  common_probes <- intersect(
    epicv2unique$Probe_base,
    v1_probes
  )
  # --- Restrict EPICv2 (use selected best probes) ---
  epicv2_overlap <- epicv2unique %>%
    filter(Probe_base %in% common_probes)
  # --- Restrict EPICv1 dataset (lazy, no collect yet) ---
  epicv1_overlap <- epicv1_vals %>%
    filter(ProbeID %in% common_probes)
  # --- Optional sanity checks ---
  message("Number of shared probes: ", length(common_probes))
  message("EPICv2 selected probes: ", nrow(epicv2_overlap))
  
  stopifnot(length(common_probes) == nrow(epicv2_overlap))
  rm(epicv1_vals, v1_probes)
  gc()
}


##to combine EPICv1 and EPICv2
if (all(c("EPICv1", "EPICv2") %in% array_types)) {
  for (data_type in data_types) {
    message("Processing: ", data_type)
    
    # --- Build file paths ---
    epicv2_vals_fn <- file.path(
      out_dir,
      paste0(out_pref, "-IlluminaHumanMethylationEPICv2-methyl-", data_type, ".parquet")
    )
    epicv1_vals_fn <- file.path(
      out_dir,
      paste0(out_pref, "-IlluminaHumanMethylationEPIC-methyl-", data_type, ".parquet")
    )
    
    # --- Open datasets ---
    epicv2_vals <- open_dataset(epicv2_vals_fn)
    epicv1_vals <- open_dataset(epicv1_vals_fn)
    
    # --- Filter datasets ---
    epicv1_filtered <- epicv1_vals %>%
      filter(ProbeID %in% common_probes)
    
    epicv2_filtered <- epicv2_vals %>%
      filter(ProbeID %in% epicv2_overlap$ProbeID)
    
    # --- Collect (now manageable size) ---
    df_v1 <- epicv1_filtered %>% collect()
    df_v2 <- epicv2_filtered %>% collect()
    
    # --- Harmonize EPICv2 IDs to Probe_base (to match v1) ---
    df_v2 <- df_v2 %>%
      mutate(ProbeID = sub("_.*$", "", ProbeID))
    # --- Align order (important!) ---
    df_v1 <- df_v1 %>%
      arrange(ProbeID)
    df_v2 <- df_v2 %>%
      arrange(ProbeID)
    # --- Optional safety check ---
    stopifnot(identical(df_v1$ProbeID, df_v2$ProbeID))
    
    # --- Combine ---
    combined_df <- bind_rows(
      df_v1,
      df_v2
    )
    
    # --- Output filename ---
    out_fn <- file.path(
      out_dir,
      paste0(out_pref, "-IlluminaHumanMethylationEPICv1-EPICv2-methyl-", data_type, ".parquet")
    )
    
    # --- Save ---
    write_parquet(combined_df, out_fn)
    message("Saved: ", out_fn)
    rm(df_v1, df_v2, combined_df, epicv2_filtered, epicv2_vals, epicv1_filtered, epicv1_vals)
    gc()
  }
}

# Add 450k restriction using existing common_probes
if (all(c("EPICv1", "EPICv2", "450k") %in% array_types)) {
  
  # --- File path ---
  k450_vals_fn <- file.path(
    out_dir,
    paste0(out_pref, "-IlluminaHumanMethylation450k-methyl-p-values.parquet")
  )
  # --- Open dataset ---
  k450_vals <- open_dataset(k450_vals_fn)
  # --- Get 450k probes ---
  k450_probes <- k450_vals %>%
    select(ProbeID) %>%
    distinct() %>%
    collect() %>%
    pull(ProbeID)
  # --- Intersect with EXISTING common_probes ---
  common_probes <- intersect(common_probes, k450_probes)
  # --- Optional diagnostics ---
  message("Shared probes across EPICv1 + EPICv2 + 450k: ", length(common_probes))
  rm(k450_vals, k450_probes)
  gc()
}

## to combine (EPICv1 + EPICv2) with 450k
if (all(c("EPICv1", "EPICv2", "450k") %in% array_types)) {
  
  for (data_type in data_types) {
    
    message("Processing: ", data_type)
    
    # --- Build file paths ---
    v1v2_vals_fn <- file.path(
      out_dir,
      paste0(out_pref, "-IlluminaHumanMethylationEPICv1-EPICv2-methyl-", data_type, ".parquet")
    )
    
    k450_vals_fn <- file.path(
      out_dir,
      paste0(out_pref, "-IlluminaHumanMethylation450k-methyl-", data_type, ".parquet")
    )
    
    # --- Open datasets ---
    v1v2_vals <- open_dataset(v1v2_vals_fn)
    k450_vals <- open_dataset(k450_vals_fn)
    
    # --- Filter datasets ---
    v1v2_filtered <- v1v2_vals %>%
      filter(ProbeID %in% common_probes)
    
    k450_filtered <- k450_vals %>%
      filter(ProbeID %in% common_probes)
    
    # --- Collect ---
    df_v1v2 <- v1v2_filtered %>% collect()
    df_450k <- k450_filtered %>% collect()
    
    # --- Align order ---
    df_v1v2 <- df_v1v2 %>%
      arrange(ProbeID)
    
    df_450k <- df_450k %>%
      arrange(ProbeID)
    
    # --- Safety check ---
    stopifnot(all(df_450k$ProbeID %in% df_v1v2$ProbeID))
    
    # Optional strict check (only if identical rows expected):
    # stopifnot(identical(df_v1v2$ProbeID, df_450k$ProbeID))
    
    # --- Combine ---
    combined_df <- bind_rows(
      df_v1v2,
      df_450k %>% mutate(Array = "450k")   # optional but recommended
    )
    
    # --- Output filename ---
    out_fn <- file.path(
      out_dir,
      paste0(out_pref, "IlluminaHumanMethylationEPICv1-EPICv2-450k-methyl-", data_type, ".parquet")
    )
    # --- Save ---
    write_parquet(combined_df, out_fn)
    message("Saved: ", out_fn)
    rm(df_v1v2, df_450k, combined_df) #remove large objects
    gc()
  }
}

if (all(c("EPICv1", "450k") %in% array_types) && !("EPICv2" %in% array_types)) 
  {
  
  # --- File paths ---
  epicv1_vals_fn <- file.path(
    out_dir,
    paste0(out_pref, "-IlluminaHumanMethylationEPIC-methyl-p-values.parquet")
  )
  
  k450_vals_fn <- file.path(
    out_dir,
    paste0(out_pref, "-IlluminaHumanMethylation450k-methyl-p-values.parquet")
  )
  
  # --- Open datasets (lazy) ---
  epicv1_vals <- open_dataset(epicv1_vals_fn)
  k450_vals   <- open_dataset(k450_vals_fn)
  
  # --- Get probe sets ---
  v1_probes <- epicv1_vals %>%
    select(ProbeID) %>%
    distinct() %>%
    collect() %>%
    pull(ProbeID)
  
  k450_probes <- k450_vals %>%
    select(ProbeID) %>%
    distinct() %>%
    collect() %>%
    pull(ProbeID)
  
  # --- Find common probes ---
  common_probes <- intersect(v1_probes, k450_probes)
  
  # --- Restrict datasets (lazy) ---
  epicv1_overlap <- epicv1_vals %>%
    filter(ProbeID %in% common_probes)
  
  k450_overlap <- k450_vals %>%
    filter(ProbeID %in% common_probes)
  
  # --- Optional sanity checks ---
  message("Number of shared probes: ", length(common_probes))
  
  # Collect counts for verification
  v1_n <- epicv1_overlap %>% select(ProbeID) %>% distinct() %>% collect() %>% nrow()
  k450_n <- k450_overlap %>% select(ProbeID) %>% distinct() %>% collect() %>% nrow()
  
  message("EPICv1 overlap probes: ", v1_n)
  message("450k overlap probes: ", k450_n)
  
  stopifnot(length(common_probes) == v1_n)
  stopifnot(length(common_probes) == k450_n)
}

## to combine EPICv1 and 450k
if (all(c("EPICv1", "450k") %in% array_types) &&
    !("EPICv2" %in% array_types)) 
  {
  for (data_type in data_types) {
    message("Processing: ", data_type)
    
    # --- Build file paths ---
    epicv1_vals_fn <- file.path(
      out_dir,
      paste0(out_pref, "-IlluminaHumanMethylationEPIC-methyl-", data_type, ".parquet")
    )
    
    k450_vals_fn <- file.path(
      out_dir,
      paste0(out_pref, "-IlluminaHumanMethylation450k-methyl-", data_type, ".parquet")
    )
    
    # --- Open datasets ---
    epicv1_vals <- open_dataset(epicv1_vals_fn)
    k450_vals   <- open_dataset(k450_vals_fn)
    
    # --- Filter datasets ---
    epicv1_filtered <- epicv1_vals %>%
      filter(ProbeID %in% common_probes)
    
    k450_filtered <- k450_vals %>%
      filter(ProbeID %in% common_probes)
    
    # --- Collect (now manageable size) ---
    df_v1 <- epicv1_filtered %>% collect()
    df_450k <- k450_filtered %>% collect()
    
    # --- Align order (important!) ---
    df_v1 <- df_v1 %>%
      arrange(ProbeID)
    
    df_450k <- df_450k %>%
      arrange(ProbeID)
    
    # --- Optional safety check ---
    stopifnot(identical(df_v1$ProbeID, df_450k$ProbeID))
    
    # --- Combine ---
    combined_df <- bind_rows(
      df_v1,
      df_450k
    )
    
    # --- Output filename ---
    out_fn <- file.path(
      out_dir,
      paste0(out_pref, "IlluminaHumanMethylationEPICv1-450k-methyl-", data_type, ".parquet")
    )
    
    # --- Save ---
    write_parquet(combined_df, out_fn)
    
    message("Saved: ", out_fn)
    rm(combined_df, df_450k, df_v1) #remove large objects 
    gc()
  }
}

# EPICv2 + 450k (no EPICv1)
if (all(c("EPICv2", "450k") %in% array_types) &&
    !("EPICv1" %in% array_types)) 
    {
  
  # --- File paths ---
  k450_vals_fn <- file.path(
    out_dir,
    paste0(out_pref, "-IlluminaHumanMethylation450k-methyl-p-values.parquet")
  )
  
  # --- Open datasets (lazy) ---
  k450_vals <- open_dataset(k450_vals_fn)
  
  # --- Get 450k probes (reference set) ---
  k450_probes <- k450_vals %>%
    select(ProbeID) %>%
    distinct() %>%
    collect() %>%
    pull(ProbeID)
  
  # --- Find common probes (via Probe_base mapping) ---
  common_probes <- intersect(
    epicv2unique$Probe_base,   # EPICv2 base IDs
    k450_probes                # 450k IDs
  )
  
  # --- Restrict EPICv2 (keep best probes per base) ---
  epicv2_overlap <- epicv2unique %>%
    filter(Probe_base %in% common_probes)
  
  # --- Restrict 450k dataset (lazy) ---
  k450_overlap <- k450_vals %>%
    filter(ProbeID %in% common_probes)
  
  # --- Optional sanity checks ---
  message("Number of shared probes: ", length(common_probes))
  message("EPICv2 selected probes: ", nrow(epicv2_overlap))
  
  stopifnot(length(common_probes) == nrow(epicv2_overlap))
}


## to combine EPICv2 and 450k
if (all(c("EPICv2", "450k") %in% array_types) &&
    !("EPICv1" %in% array_types)) 
  {
  
  for (data_type in data_types) {
    
    message("Processing: ", data_type)
    
    # --- Build file paths ---
    epicv2_vals_fn <- file.path(
      out_dir,
      paste0(out_pref, "-IlluminaHumanMethylationEPICv2-methyl-", data_type, ".parquet")
    )
    
    k450_vals_fn <- file.path(
      out_dir,
      paste0(out_pref, "-IlluminaHumanMethylation450k-methyl-", data_type, ".parquet")
    )
    
    # --- Open datasets ---
    epicv2_vals <- open_dataset(epicv2_vals_fn)
    k450_vals   <- open_dataset(k450_vals_fn)
    
    # --- Filter datasets ---
    k450_filtered <- k450_vals %>%
      filter(ProbeID %in% common_probes)
    
    epicv2_filtered <- epicv2_vals %>%
      filter(ProbeID %in% epicv2_overlap$ProbeID)
    
    # --- Collect ---
    df_450k <- k450_filtered %>% collect()
    df_v2   <- epicv2_filtered %>% collect()
    
    # --- Harmonize EPICv2 IDs to base CpG IDs ---
    df_v2 <- df_v2 %>%
      mutate(ProbeID = sub("_.*$", "", ProbeID))
    
    # --- Align order ---
    df_450k <- df_450k %>%
      arrange(ProbeID)
    
    df_v2 <- df_v2 %>%
      arrange(ProbeID)
    
    # --- Safety check ---
    stopifnot(identical(df_450k$ProbeID, df_v2$ProbeID))
    
    # --- Combine ---
    combined_df <- bind_rows(
      df_450k,
      df_v2
    )
    
    # --- Output filename ---
    out_fn <- file.path(
      out_dir,
      paste0(out_pref, "IlluminaHumanMethylationEPICv2-450k-methyl-", data_type, ".parquet")
    )
    
    # --- Save ---
    write_parquet(combined_df, out_fn)
    
    message("Saved: ", out_fn)
    rm(combined_df, df_v2, df_450k)
    gc()
  }
}


