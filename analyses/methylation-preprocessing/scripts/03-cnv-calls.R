# use minfi mset object to call CNVs using conumee2
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
  make_option(opt_str = "--base_dir", type = "character", default = NULL,
              help = "The absolute path of the base directory containing sample 
              array IDAT files.",
              metavar = "character"),
  make_option(opt_str = "--output_basename", type = "character", default = NULL,
              help = "The same output basename for 01-preprocess-illumina-arrays",
              metavar = "character"),
  make_option(opt_str = "--manifest_file", type = "character",
              help = "Input manifest file with 'file_name' and
              'Bioassay_ID' columns"),
  make_option(opt_str = "--n_cores", type = 'integer',
              default=1, help="number of cores for parallelisation of minfi::detectionP. Default is 1")
)


# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
base_dir <- opt$base_dir
manifest_file <- opt$manifest_file
n_cores <- opt$n_cores
out_base <- opt$output_basename


base_dir <- 'sorted_idats_output_dir/IlluminaHumanMethylationEPIC'
n_cores <- 4
manifest_file <- 'controls_and_dicer_manifest.tsv'
out_base <- 'test-out/'
dataset <- basename(base_dir)


m_set_file <- paste0(out_base, "-", dataset, '-m-set.qs2')
GRset <- qs_read(m_set_file)


anno <- CNV.create_anno(array_type = "EPIC") 

# read manifest to obtain the IDAT prefix from the `file_name` and its matched `Bioassay_ID` column
man_df <- read_tsv(file = manifest_file, show_col_types = FALSE) %>% 
  dplyr::select(file_name, Bioassay_ID, sample_type, platform) %>%
  dplyr::mutate(file_name = gsub("(_Red|_Grn).*", "", file_name)) %>%
  dplyr::mutate(file_name = basename(file_name)) %>%
  dplyr::filter(platform == "HumanMethylationEPIC") %>%
  unique()

normals <- man_df %>% 
  dplyr::filter(sample_type == 'Normal') %>%
  dplyr::select(file_name)

reference_samples <- intersect(normals$file_name, colnames(GRset))



query <- CNV.load(
  GRset[, '204305120005_R08C01'],
  anno = anno
)



cnv <- conumee2::CNV.fit(
  GRset[, '204305120005_R08C01'],
  GRset[, reference_samples],
  anno
)
