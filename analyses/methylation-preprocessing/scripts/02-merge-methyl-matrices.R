# Merge per-array methylation matrices by ProbeID, retaining samples from
# every input batch and intersecting probes across the selected arrays.
# Jessica Daggett, 05/12/2026

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(qs2))

option_list <- list(
  make_option("--output_dir", type = "character", default = NULL,
              help = "Directory containing input and output parquet files"),
  make_option("--output_prefix", type = "character", default = NULL,
              help = "Output file prefix"),
  make_option("--input_prefix", type = "character", default = NULL,
              help = "Comma-separated input file prefixes; defaults to output_prefix")
)
opt <- parse_args(OptionParser(option_list = option_list))
out_dir <- opt$output_dir
out_pref <- opt$output_prefix
if (is.null(out_dir) || !dir.exists(out_dir)) {
  stop("--output_dir must name an existing directory.")
}
if (is.null(out_pref) || !nzchar(trimws(out_pref))) {
  stop("--output_prefix is required.")
}
input_prefixes <- if (is.null(opt$input_prefix)) {
  out_pref
} else {
  unique(trimws(strsplit(opt$input_prefix, ",", fixed = TRUE)[[1]]))
}
if (length(input_prefixes) == 0L || any(!nzchar(input_prefixes))) {
  stop("--input_prefix must contain nonempty prefixes.")
}

data_types <- c("beta-values-masked", "m-values-masked",
                "m-values-unmasked", "cn-values")
input_data_types <- c("p-values", data_types)
array_names <- c(
  EPICv1 = "IlluminaHumanMethylationEPIC",
  EPICv2 = "IlluminaHumanMethylationEPICv2",
  `450k` = "IlluminaHumanMethylation450k"
)

# Match exact per-array filenames. Prefix text and old merged outputs cannot
# introduce an array type or accidentally select an additional batch.
input_files <- lapply(array_names, function(array_name) {
  batches <- lapply(input_prefixes, function(prefix) {
    paths <- setNames(file.path(out_dir, paste0(
      prefix, "-", array_name, "-methyl-", input_data_types, ".parquet"
    )), input_data_types)
    present <- file.exists(paths)
    if (!any(present)) return(NULL)
    if (!all(present)) {
      stop("Incomplete inputs for ", prefix, " / ", array_name,
           ". Missing: ", paste(basename(paths[!present]), collapse = ", "))
    }
    paths
  })
  batches <- Filter(Negate(is.null), batches)
  setNames(lapply(input_data_types, function(data_type) {
    vapply(batches, function(paths) paths[[data_type]], character(1))
  }), input_data_types)
})
array_types <- names(input_files)[vapply(
  input_files, function(x) length(x[["p-values"]]) > 0L, logical(1)
)]
if (length(array_types) == 0L) {
  stop("No per-array methylation outputs found for prefixes '",
       paste(input_prefixes, collapse = ", "), "' in ", out_dir, ".")
}
message("Array types: ", paste(array_types, collapse = ", "))

check_probe_ids <- function(ids, path) {
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("ProbeID must be nonmissing, nonempty and unique in ", basename(path))
  }
}
read_probe_ids <- function(path) {
  ids <- open_dataset(path) %>% select(ProbeID) %>% collect() %>% pull(ProbeID)
  check_probe_ids(ids, path)
  ids
}

# Validate sample schemas before writing outputs. Reused sample names must be
# resolved explicitly rather than silently acquiring .x/.y join suffixes.
seen_samples <- character()
array_probes <- list()
for (array_type in array_types) {
  paths <- input_files[[array_type]]
  probe_sets <- list()
  for (batch in seq_along(paths[["p-values"]])) {
    batch_samples <- NULL
    for (data_type in input_data_types) {
      path <- paths[[data_type]][batch]
      columns <- colnames(open_dataset(path))
      if (!"ProbeID" %in% columns || anyDuplicated(columns)) {
        stop("Expected unique columns including ProbeID in ", basename(path))
      }
      samples <- setdiff(columns, "ProbeID")
      if (length(samples) == 0L) stop("No sample columns in ", basename(path))
      if (is.null(batch_samples)) {
        batch_samples <- samples
      } else if (!setequal(samples, batch_samples)) {
        stop("Sample columns differ between measurement types in ", basename(path))
      }
    }
    overlap <- intersect(seen_samples, batch_samples)
    if (length(overlap) > 0L) {
      stop("Sample names occur in multiple input batches/arrays: ",
           paste(overlap, collapse = ", "))
    }
    seen_samples <- c(seen_samples, batch_samples)
    probe_sets[[batch]] <- read_probe_ids(paths[["p-values"]][batch])
  }
  common <- Reduce(intersect, probe_sets)
  if (length(common) == 0L) stop("No common probes across ", array_type, " batches.")
  for (batch in seq_along(probe_sets)) {
    message(array_type, " batch ", batch, ": retaining ", length(common),
            " probes; excluding ", length(probe_sets[[batch]]) - length(common), ".")
  }
  array_probes[[array_type]] <- common
}

if (length(array_types) == 1L &&
    length(input_files[[array_types]][["p-values"]]) == 1L) {
  message("Only one array batch found; leaving per-array outputs unchanged.")
  quit(save = "no", status = 0)
}

# For cross-array merges, select one EPICv2 probe per base CpG using all
# EPICv2 samples. Only collect duplicate sites, and combine sums/counts so
# unequal batch sizes and missing p-values are weighted correctly.
epicv2_mapping <- NULL
if (length(array_types) > 1L && "EPICv2" %in% array_types) {
  ids <- array_probes[["EPICv2"]]
  bases <- sub("_.*$", "", ids)
  duplicate_bases <- unique(bases[duplicated(bases)])
  duplicate_ids <- ids[bases %in% duplicate_bases]
  p_sums <- numeric(length(duplicate_ids))
  p_counts <- numeric(length(duplicate_ids))
  if (length(duplicate_ids) > 0L) {
    for (path in input_files[["EPICv2"]][["p-values"]]) {
      detP <- open_dataset(path) %>%
        filter(ProbeID %in% duplicate_ids) %>% collect()
      detP <- detP[match(duplicate_ids, detP$ProbeID), ]
      values <- as.matrix(detP[, setdiff(names(detP), "ProbeID"), drop = FALSE])
      p_sums <- p_sums + rowSums(values, na.rm = TRUE)
      p_counts <- p_counts + rowSums(!is.na(values))
      rm(detP, values)
    }
  }
  scores <- rep(Inf, length(duplicate_ids))
  has_values <- p_counts > 0
  scores[has_values] <- p_sums[has_values] / p_counts[has_values]
  best <- tibble(ProbeID = duplicate_ids,
                 Probe_base = sub("_.*$", "", duplicate_ids),
                 detP_mean = scores) %>%
    arrange(Probe_base, detP_mean, ProbeID) %>%
    group_by(Probe_base) %>% slice_head(n = 1) %>% ungroup() %>%
    select(ProbeID, Probe_base)
  epicv2_mapping <- bind_rows(
    best,
    tibble(ProbeID = ids[!bases %in% duplicate_bases],
           Probe_base = bases[!bases %in% duplicate_bases])
  ) %>% arrange(Probe_base)
  qs_save(epicv2_mapping,
          file.path(out_dir, paste0(out_pref, "-EPICv2-unique-mapping.qs2")))
  array_probes[["EPICv2"]] <- epicv2_mapping$Probe_base
  gc()
}

# Load each batch separately: an Arrow dataset containing multiple files
# stacks rows, whereas these wide matrices need samples joined as columns.
load_array_values <- function(array_type, data_type, common_probes) {
  selected_ids <- common_probes
  if (array_type == "EPICv2" && !is.null(epicv2_mapping)) {
    selected_ids <- epicv2_mapping$ProbeID[
      match(common_probes, epicv2_mapping$Probe_base)
    ]
  }
  combined <- NULL
  for (path in input_files[[array_type]][[data_type]]) {
    values <- open_dataset(path) %>%
      filter(ProbeID %in% selected_ids) %>% collect()
    check_probe_ids(values$ProbeID, path)
    if (!setequal(values$ProbeID, selected_ids)) {
      stop("Missing selected probes in ", basename(path))
    }
    if (array_type == "EPICv2" && !is.null(epicv2_mapping)) {
      values$ProbeID <- sub("_.*$", "", values$ProbeID)
    }
    values <- values %>% arrange(ProbeID)
    combined <- if (is.null(combined)) values else {
      inner_join(combined, values, by = "ProbeID", relationship = "one-to-one")
    }
    rm(values)
  }
  combined
}

merge_arrays <- function(selected_arrays) {
  # Sort by radix: dplyr's arrange() orders character columns in the C locale,
  # and the identical() check below compares against a frame ordered that way.
  # Plain sort() would follow LC_COLLATE and could disagree on the '.' and
  # uppercase characters in ch./rs probe names under some locales.
  common_probes <- sort(Reduce(intersect, array_probes[selected_arrays]),
                        method = "radix")
  if (length(common_probes) == 0L) {
    stop("No common probes across ", paste(selected_arrays, collapse = ", "), ".")
  }
  message("Combining ", paste(selected_arrays, collapse = ", "),
          ": ", length(common_probes), " shared probes.")
  for (array_type in selected_arrays) {
    message(array_type, ": excluding ",
            length(array_probes[[array_type]]) - length(common_probes),
            " probes outside the shared set.")
  }
  # Same-array batches keep their original IDs and use an explicit merged
  # suffix so the result cannot be rediscovered as a per-array input on reruns.
  output_array <- if (length(selected_arrays) == 1L) {
    paste0(array_names[[selected_arrays]], "-merged")
  } else {
    paste0("IlluminaHumanMethylation", paste(selected_arrays, collapse = "-"))
  }
  expected_samples <- unlist(lapply(selected_arrays, function(array_type) {
    unlist(lapply(input_files[[array_type]][["p-values"]], function(path) {
      setdiff(colnames(open_dataset(path)), "ProbeID")
    }), use.names = FALSE)
  }), use.names = FALSE)
  for (data_type in data_types) {
    message("Processing: ", data_type)
    combined <- NULL
    for (array_type in selected_arrays) {
      values <- load_array_values(array_type, data_type, common_probes)
      combined <- if (is.null(combined)) values else {
        inner_join(combined, values, by = "ProbeID", relationship = "one-to-one")
      }
      rm(values)
    }
    stopifnot(identical(combined$ProbeID, common_probes),
              setequal(setdiff(names(combined), "ProbeID"), expected_samples),
              ncol(combined) == length(expected_samples) + 1L)
    out_fn <- file.path(out_dir, paste0(
      out_pref, "-", output_array, "-methyl-", data_type, ".parquet"
    ))
    write_parquet(combined, out_fn, compression = "zstd")
    message("Saved: ", out_fn)
    rm(combined)
    gc()
  }
}

# Preserve the EPICv1+EPICv2 output as well as the three-array output when
# all arrays are present. All combinations share the same implementation.
if (all(c("EPICv1", "EPICv2", "450k") %in% array_types)) {
  merge_arrays(c("EPICv1", "EPICv2"))
}
merge_arrays(array_types)
