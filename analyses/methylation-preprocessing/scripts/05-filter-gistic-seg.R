#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage: 05-filter-gistic-seg.R <input.seg> <bioassay_ids.txt> <output.seg>\n",
    "\n",
    "Keeps GISTIC SEG rows whose Sample value is listed in bioassay_ids.txt.\n",
    sep = ""
  )
}

if (length(args) == 1L && args[[1L]] %in% c("-h", "--help")) {
  usage()
  quit(status = 0L)
}

if (length(args) != 3L) {
  usage()
  quit(status = 2L)
}

input_seg <- args[[1L]]
bioassay_id_file <- args[[2L]]
output_seg <- args[[3L]]

for (path in c(input_seg, bioassay_id_file)) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }
}

bioassay_ids <- unique(trimws(readLines(bioassay_id_file, warn = FALSE)))
bioassay_ids <- bioassay_ids[nzchar(bioassay_ids)]
if (!length(bioassay_ids)) {
  stop("No bioassay IDs found in: ", bioassay_id_file, call. = FALSE)
}

seg <- read.delim(
  input_seg,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!"Sample" %in% names(seg)) {
  stop("SEG file must contain a Sample column: ", input_seg, call. = FALSE)
}

filtered_seg <- seg[seg$Sample %in% bioassay_ids, , drop = FALSE]
matched_ids <- unique(filtered_seg$Sample)

if (!nrow(filtered_seg)) {
  stop("None of the requested bioassay IDs are present in: ", input_seg, call. = FALSE)
}

output_dir <- dirname(output_seg)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
write.table(
  filtered_seg,
  output_seg,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "Wrote ", nrow(filtered_seg), " segments for ", length(matched_ids),
  " of ", length(bioassay_ids), " requested bioassay IDs to ", output_seg
)
