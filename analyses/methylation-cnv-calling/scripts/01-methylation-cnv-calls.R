# Use a Noob-normalized MethylSet to call CNVs using conumee2.
# Jessica Daggett
# 05/07/2026

Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS = 1
)

# Load libraries:
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(arrow))
suppressWarnings(
  suppressPackageStartupMessages(library(minfi))
)

suppressPackageStartupMessages(library("conumee2"))


# Magrittr pipe
`%>%` <- dplyr::`%>%`

# set up optparse options
option_list <- list(
  make_option(opt_str = "--rg_set_file", type = "character", default = NULL,
              help = "QC-filtered RGChannelSet serialized as a .qs2 file.",
              metavar = "file"),
  make_option(opt_str = "--output_basename", type = "character", default = NULL,
              help = "Path and prefix for the output segment files.",
              metavar = "character"),
  make_option(opt_str = "--manifest_file", type = "character",
              help = "Input manifest file with 'file_name' and
              'Bioassay_ID' columns"),
  make_option(opt_str = "--array_type", type = 'character',
              default="EPIC", help="Array type: EPIC, EPICv2, or 450k.")
)


# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
rg_set_file <- opt$rg_set_file
manifest_file <- opt$manifest_file
out_base <- opt$output_basename
array_type <- opt$array_type

if (is.null(rg_set_file) || !file.exists(rg_set_file)) {
  stop("--rg_set_file must name an existing .qs2 RGset file.")
}
if (is.null(manifest_file) || !file.exists(manifest_file)) {
  stop("--manifest_file must name an existing manifest file.")
}
if (is.null(out_base) || !nzchar(out_base)) {
  stop("--output_basename is required.")
}

# Normalize known platform naming variants before filtering. This accepts the
# current manifest spellings (for example, HumanMethylationEPICv2 and Illumina
# Infinium HumanMethylationEPICv2) as well as harmless spacing/punctuation and
# v1 version suffixes.
platform_matches_array <- function(platform, array_type) {
  platform_key <- tolower(gsub("[^[:alnum:]]", "", platform))

  switch(
    array_type,
    EPICv2 = grepl("humanmethylationepicv?2", platform_key),
    EPIC = grepl("humanmethylationepic($|v?1)", platform_key),
    `450k` = grepl("humanmethylation(450k|450)", platform_key),
    HM450 = grepl("humanmethylation(450k|450)", platform_key),
    `450` = grepl("humanmethylation(450k|450)", platform_key),
    stop("Unsupported --array_type: ", array_type,
         ". Expected EPIC, EPICv2, or 450k.", call. = FALSE)
  )
}

# Read the manifest and retain the matching array platform.
manifest_df <- read_tsv(file = manifest_file, show_col_types = FALSE) %>%
  dplyr::select(file_name, Bioassay_ID, sample_type, platform) %>%
  dplyr::mutate(file_name = gsub("(_Red|_Grn).*", "", file_name)) %>%
  dplyr::mutate(file_name = basename(file_name))

man_df <- manifest_df %>%
  dplyr::filter(platform_matches_array(platform, array_type)) %>%
  unique()


normals <- man_df %>% 
  dplyr::filter(sample_type == 'Normal') %>%
  dplyr::select(Bioassay_ID)

if (nrow(man_df) == 0) {
  available_platforms <- sort(unique(manifest_df$platform))
  stop(
    "No manifest rows matched array type '", array_type, "'. Observed platform values: ",
    paste(available_platforms, collapse = ", "),
    call. = FALSE
  )
}

# Load the QC-filtered RGChannelSet and apply Noob normalization here.  Conumee2
# segmentation needs the intensity signal that funnorm/quantile normalization
# does not preserve.
RGset <- qs_read(rg_set_file)
MSet <- minfi::preprocessNoob(RGset)
rm(RGset)
gc()

# Make MSet names Bioassay IDs.
MSet <- MSet[, colnames(MSet) %in% man_df$file_name]
name_map <- setNames(man_df$Bioassay_ID, man_df$file_name)
colnames(MSet) <- name_map[colnames(MSet)]

##get reference vs query samples 
sample_names <- colnames(MSet)
reference_samples <- intersect(normals$Bioassay_ID, sample_names)
query_samples     <- setdiff(sample_names, reference_samples)
if (length(reference_samples) == 0) {
  stop("No normal samples in the RGset matched the supplied manifest.")
}
if (length(query_samples) == 0) {
  stop("No non-normal samples in the RGset matched the supplied manifest.")
}
#get ref vs query mset
MSet_ref   <- MSet[, reference_samples]
MSet_query <- MSet[, query_samples]
#load both cnvs
ref   <- CNV.load(MSet_ref)
query <- CNV.load(MSet_query)


## get annotation - epicv2 uses hg38, others use hg19. Perform liftover later
if (array_type == "EPICv2") {
  message("Creating annotation for EPICv2 (hg38)")
  anno <- CNV.create_anno(
    array_type = "EPICv2",
    genome = "hg38"
  )
} else if (array_type == "EPIC") {
  
  message("Creating annotation for EPIC (hg19)")
  
  anno <- CNV.create_anno(
    array_type = "EPIC"
  )
} else if (array_type %in% c("450k", "HM450", "450")) {
  message("Creating annotation for 450k (hg19)")
  anno <- CNV.create_anno(
    array_type = "450k"
  )
} else {
  warning("Unknown array type — defaulting to EPIC hg19")
  anno <- CNV.create_anno(
    array_type = "EPIC"
  )
  
}
#fit cnvs
x <- CNV.fit(query, ref, anno)
x <- CNV.bin(x)
#x <- CNV.detail(x) #only need if you provide detail regions 

# Get segments; parameters are tuned for array type.
if (array_type %in% c("EPIC", "EPICv2")) {
  
  message("Using EPIC-optimized segmentation parameters")
  
  x <- CNV.segment(
    x,
    alpha = 0.001,
    nperm = 50000,
    min.width = 5,
    undo.splits = "sdundo",
    undo.SD = 2.5
  )
  
} else {
  
  message("Using 450k default segmentation parameters")
  
  x <- CNV.segment(
    x,
    alpha = 0.001,
    nperm = 50000,
    min.width = 5,
    undo.splits = "sdundo",
    undo.SD = 2.2
  )
}


segments_file <- paste0(out_base, "-segments.seg")
gistic_file <- paste0(out_base, "-gistic.seg")


## these will form the input files to gistic, after some post processing to get the headers in the right format
segments <- CNV.write(x, what = "segments")
segments_df <- bind_rows(segments)

data.table::fwrite(
  segments_df,
  file = segments_file,
  sep = "\t",
  col.names = TRUE
)

# Write GISTIC input.
CNV.write(x, what="gistic", file=gistic_file)

