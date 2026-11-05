# QC of the GISTIC cohort.  Run before 04-run-gistic.sh: reports, per array,
# how many cohort samples contribute, how amplitude and segmentation compare
# across arrays, and how many segments survive liftover to hg38.  The cohort is
# the samples in --bioassay_ids, or all samples when it is not given.
library(optparse)

option_list <- list(
  make_option(opt_str = "--seg_dir", type = "character", default = NULL,
              help = "Directory with per-array *-gistic.seg files and combined_hg38.gistic.seg",
              metavar = "character"),
  make_option(opt_str = "--threshold", type = "double", default = 0.2,
              help = "|Seg.CN| threshold for an altered segment; matches GISTIC -ta/-td [default %default]",
              metavar = "double"),
  make_option(opt_str = "--bioassay_ids", type = "character", default = NULL,
              help = "Optional file of Bioassay_IDs (one per line). Summarize only these samples.",
              metavar = "character"),
  make_option(opt_str = "--output_file", type = "character", default = NULL,
              help = "Output TSV. Default: <seg_dir>/gistic-cohort-qc.tsv",
              metavar = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$seg_dir) || !dir.exists(opt$seg_dir)) {
  stop("--seg_dir must be an existing directory.")
}
output_file <- if (is.null(opt$output_file)) {
  file.path(opt$seg_dir, "gistic-cohort-qc.tsv")
} else {
  opt$output_file
}

cohort_ids <- NULL
if (!is.null(opt$bioassay_ids)) {
  if (!file.exists(opt$bioassay_ids)) stop("Bioassay ID file not found: ", opt$bioassay_ids)
  cohort_ids <- unique(trimws(readLines(opt$bioassay_ids, warn = FALSE)))
  cohort_ids <- cohort_ids[nzchar(cohort_ids)]
  if (!length(cohort_ids)) stop("No bioassay IDs found in: ", opt$bioassay_ids)
}

combined_file <- file.path(opt$seg_dir, "combined_hg38.gistic.seg")
if (!file.exists(combined_file)) {
  stop("combined_hg38.gistic.seg not found in: ", opt$seg_dir)
}
combined <- read.delim(combined_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!is.null(cohort_ids)) {
  combined <- combined[combined$Sample %in% cohort_ids, , drop = FALSE]
}

seg_files <- list.files(opt$seg_dir, pattern = "gistic\\.seg$", full.names = TRUE)
seg_files <- seg_files[basename(seg_files) != "combined_hg38.gistic.seg"]
if (!length(seg_files)) {
  stop("No per-array GISTIC SEG files found in: ", opt$seg_dir)
}

# Same filename convention as 02-liftover-methylation-segments.R
detect_platform <- function(f) {
  f <- basename(f)
  if (grepl("EPICv2", f, ignore.case = TRUE)) return("EPICv2")
  if (grepl("EPIC-", f, ignore.case = TRUE)) return("EPIC")
  if (grepl("450", f, ignore.case = TRUE)) return("450k")
  "unknown"
}

summarize_array <- function(f) {
  seg <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(cohort_ids)) {
    seg <- seg[seg$Sample %in% cohort_ids, , drop = FALSE]
  }
  if (!nrow(seg)) {
    return(NULL)
  }
  width <- seg$End_Position - seg$Start_Position + 1
  altered <- abs(seg$Seg.CN) >= opt$threshold
  samples <- unique(seg$Sample)
  n_combined <- sum(combined$Sample %in% samples)
  q <- quantile(seg$Seg.CN, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)
  per_sample <- table(seg$Sample)

  data.frame(
    platform = detect_platform(f),
    source_file = basename(f),
    n_samples = length(samples),
    n_segments = nrow(seg),
    median_segments_per_sample = median(as.numeric(per_sample)),
    median_segment_width_mb = round(median(width) / 1e6, 2),
    median_num_markers = median(seg$Num_Markers),
    seg_cn_p05 = round(q[1], 3),
    seg_cn_p25 = round(q[2], 3),
    seg_cn_median = round(q[3], 3),
    seg_cn_p75 = round(q[4], 3),
    seg_cn_p95 = round(q[5], 3),
    seg_cn_sd = round(sd(seg$Seg.CN), 3),
    frac_segments_altered = round(mean(altered), 3),
    frac_genome_altered = round(sum(width[altered]) / sum(width), 3),
    n_segments_in_combined = n_combined,
    frac_retained_after_liftover = round(n_combined / nrow(seg), 3),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(seg_files, summarize_array))
if (is.null(summary_df)) {
  stop("None of the requested bioassay IDs are present in the per-array SEG files.")
}
if (!is.null(cohort_ids)) {
  message(sum(summary_df$n_samples), " of ", length(cohort_ids),
          " requested bioassay IDs found across arrays.")
}
summary_df <- summary_df[order(summary_df$platform), ]
if (any(summary_df$platform == "unknown")) {
  stop("Could not determine array platform for: ",
       paste(summary_df$source_file[summary_df$platform == "unknown"], collapse = ", "))
}

dup_samples <- unlist(lapply(seg_files, function(f) {
  ids <- unique(read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)$Sample)
  if (is.null(cohort_ids)) ids else ids[ids %in% cohort_ids]
}))
if (anyDuplicated(dup_samples)) {
  warning("Sample IDs appear in more than one array file: ",
          paste(unique(dup_samples[duplicated(dup_samples)]), collapse = ", "),
          call. = FALSE)
}

if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}
write.table(summary_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote segment QC summary for ", nrow(summary_df), " array file(s) to ", output_file)
