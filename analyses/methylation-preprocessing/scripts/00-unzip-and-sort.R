# Unzips and sorts idats files into array type folders
# avoids an error in minfi trying to process mixed arrays
# Largely copied minfi.readmetharray.
# Major steps:
#   Find idat files in input directory
#   unzip them
#   guess array type
#   sort into array type folders

# Alex Sickler for Pediatric OpenTargets
# 04/03/2025

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(illuminaio))
suppressPackageStartupMessages(library(BiocParallel))
suppressWarnings(
    suppressPackageStartupMessages(library(minfi))
)

.guessArrayTypes <- function(nProbes) {
    if (nProbes >= 622000 && nProbes <= 623000) {
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylation450k"
        )
    } else if (nProbes >= 1050000 && nProbes <= 1053000) {
        # NOTE: "Current EPIC scan type"
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylationEPIC"
        )
    } else if (nProbes >= 1032000 && nProbes <= 1033000) {
        # NOTE: "Old EPIC scan type"
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylationEPIC"
        )
    } else if (nProbes >= 1105000 && nProbes <= 1105300) {
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylationEPICv2"
        )
    } else if (nProbes >= 55200 && nProbes <= 55400) {
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylation27k"
        )
    } else if (nProbes >= 54700 && nProbes <= 54800) {
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylationAllergy"
        )
    } else if (nProbes >= 41000 && nProbes <= 41100) {
        arrayAnnotation <- c(
            array = "HorvathMammalMethylChip40"
        )
    } else if (nProbes >= 43650 && nProbes <= 43680) {
        arrayAnnotation <- c(
            array = "IlluminaHumanMethylationAllergy"
        )
    } else {
        arrayAnnotation <- c(array = "Unknown")
    }
    arrayAnnotation
}

# set up optparse options
option_list <- list(
    make_option(
        opt_str = "--base_dir",
        type = "character", default = NULL,
        help = "The absolute path of the base directory containing sample array IDAT files.",
        metavar = "character"
    ),
    make_option(
        opt_str = "--output_basename",
        type = "character", default = NULL,
        help = "The absolute path of the base directory containing sample array IDAT files.",
        metavar = "character"
    ),
    make_option(
        opt_str = "--manifest_file", type = "character",
        help = "Input manifest file with 'file_name' and
              'Bioassay_ID' columns. An optional 'unzipped_file_name' column
              is used when the corresponding compressed file is unavailable."
    )
)

# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
base_dir <- opt$base_dir
out_base <- opt$output_basename

manifest_df <- read_tsv(file = opt$manifest_file, show_col_types = FALSE)
required_columns <- c("file_name", "Bioassay_ID")
missing_columns <- setdiff(required_columns, names(manifest_df))
if (length(missing_columns) > 0) {
    stop("Manifest is missing required column(s): ",
         paste(missing_columns, collapse = ", "))
}

# A rerun may occur after gunzip has removed the original .gz files. Use an
# explicitly provided unzipped filename when available; otherwise derive it
# from file_name by removing a trailing .gz suffix.
unzipped_file_names <- sub("\\.gz$", "", manifest_df$file_name)
if ("unzipped_file_name" %in% names(manifest_df)) {
    supplied_unzipped_names <- manifest_df$unzipped_file_name
    has_supplied_name <- !is.na(supplied_unzipped_names) &
        nzchar(supplied_unzipped_names)
    unzipped_file_names[has_supplied_name] <-
        supplied_unzipped_names[has_supplied_name]
}

man_df <- manifest_df %>%
    transmute(file_name, unzipped_file_name = unzipped_file_names, Bioassay_ID) %>%
    unique()

message("Finding IDAT files in ", base_dir)

idat_files <- list.files(
    path = base_dir,
    pattern = "idat",
    full.names = TRUE,
    recursive = TRUE
)

out_dir <- paste0(out_base, "_output_dir")
dir.create(out_dir)

# Match the compressed filename first. If it is unavailable, use the
# uncompressed filename so the script can be rerun after the initial gunzip.
manifest_file_names <- basename(man_df$file_name)
unzipped_file_names <- basename(man_df$unzipped_file_name)
idat_file_names <- basename(idat_files)
matched_file_indices <- match(manifest_file_names, idat_file_names)
use_unzipped_file <- is.na(matched_file_indices)
matched_file_indices[use_unzipped_file] <- match(
    unzipped_file_names[use_unzipped_file], idat_file_names
)

not_in_folder <- manifest_file_names[is.na(matched_file_indices)]

out_file = paste0(out_base, "_additional_files.txt")
writeLines(
    setdiff(idat_file_names, idat_file_names[matched_file_indices]),
    file.path(out_file)
)

if (length(not_in_folder) > 0) {
    message("Error: Can't find the following files in base_dir:")
    message(paste(not_in_folder, sep = "\n"))
    stop()
}

# Retain the file selected for each manifest entry. This prefers compressed
# inputs, but falls back to their uncompressed counterparts on reruns.
idat_files <- idat_files[matched_file_indices]

message("Unzipping IDAT files")

BPREDO <- list()
BPPARAM <- SerialParam()

idat_files <- bplapply(idat_files, function(xx) {
    if (grepl(".gz$", xx)) {
        message("Unzipping ", xx)
        gunzip(xx, overwrite = TRUE)
        xx <- gsub(".gz$", "", xx)
    }
    xx
}, BPREDO = BPREDO, BPPARAM = BPPARAM)

message(idat_files)

message("Reading IDAT files")
all_quants <- bplapply(idat_files, function(xx) {
    message("Reading ", xx)
    quants <- readIDAT(xx)[["Quants"]]
}, BPREDO = BPREDO, BPPARAM = BPPARAM)

message("Sorting IDAT files into array type folders")
all_n_probes <- vapply(all_quants, nrow, integer(1L))
array_types <- cbind(do.call(rbind, lapply(all_n_probes, .guessArrayTypes)),
    size = all_n_probes
)

unique_array_types <- unique(array_types[, "array"])
for (array_type in unique_array_types) {
    message("Sorting ", array_type, " IDAT files")
    array_type_files <- idat_files[array_types[, "array"] == array_type]
    array_type_dir <- file.path(out_dir, array_type)
    if (!dir.exists(array_type_dir)) {
        dir.create(array_type_dir)
    }
    for (file in array_type_files) {
        file.copy(file.path(file), array_type_dir)
    }
}
