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
option_list <- list(make_option(
    opt_str = "--base_dir",
    type = "character", default = NULL,
    help = "The absolute path of the base directory containing sample array IDAT files.",
    metavar = "character"
))


# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
base_dir <- opt$base_dir

message("Finding IDAT files in ", base_dir)

idat_files <- list.files(
    path = base_dir,
    pattern = "idat",
    full.names = TRUE,
    recursive = TRUE
)


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

dir.create("output_dir")

unique_array_types <- unique(array_types[, "array"])
for (array_type in unique_array_types) {
    message("Sorting ", array_type, " IDAT files")
    array_type_files <- idat_files[array_types[, "array"] == array_type]
    array_type_dir <- file.path("output_dir", array_type)
    if (!dir.exists(array_type_dir)) {
        dir.create(array_type_dir)
    }
    for (file in array_type_files) {
        file.copy(file.path(file), array_type_dir)
    }
}
