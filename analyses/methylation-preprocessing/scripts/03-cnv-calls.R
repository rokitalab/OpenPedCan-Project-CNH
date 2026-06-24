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
              default=1, help="number of cores for parallelisation of minfi::detectionP. Default is 1"),
  make_option(opt_str = "--array_type", type = 'character',
              default="EPIC", help="short form array type, either EPIC, EPICv2 or 450k")
)


# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
base_dir <- opt$base_dir
manifest_file <- opt$manifest_file
n_cores <- opt$n_cores
out_base <- opt$output_basename
array_type <- opt$array_type
#for local testing
#base_dir <- 'sorted_idats_output_dir/IlluminaHumanMethylationEPICv2'
#n_cores <- 4
#manifest_file <- 'controls_and_dicer_manifest.tsv'
#out_base <- 'test-out/test'
#array_type <- 'EPICv2'

dataset <- basename(base_dir)

# read manifest to obtain the IDAT prefix from the `file_name` and its matched `Bioassay_ID` column
man_df <- read_tsv(file = manifest_file, show_col_types = FALSE) %>% 
  dplyr::select(file_name, Bioassay_ID, sample_type, platform) %>%
  dplyr::mutate(file_name = gsub("(_Red|_Grn).*", "", file_name)) %>%
  dplyr::mutate(file_name = basename(file_name)) %>%
  dplyr::filter(platform %in% c(paste0("Illumina Infinium HumanMethylation",array_type),paste0("HumanMethylation",array_type))) %>%
  unique()


normals <- man_df %>% 
  dplyr::filter(sample_type == 'Normal') %>%
  dplyr::select(Bioassay_ID)

##m set file has all data for 
m_set_file <- paste0(out_base, "-", dataset, '-m-set.qs2')
MSet <- qs_read(m_set_file)

#make mset names bioassay ids
intersect_samples <- intersect(colnames(MSet), man_df$file_name)
MSet <- MSet[, colnames(MSet) %in% man_df$file_name]
name_map <- setNames(man_df$Bioassay_ID, man_df$file_name)
colnames(MSet) <- name_map[colnames(MSet)]

##get reference vs query samples 
sample_names <- colnames(MSet)
reference_samples <- intersect(normals$Bioassay_ID, sample_names)
query_samples     <- setdiff(sample_names, reference_samples)
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

#get segments - parameters tuned for array type
x <- CNV.bin(x)

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


segments_file <- paste0(out_base, "-", dataset, "-segments",".seg")
gistic_file <- paste0(out_base, "-", dataset,"-gistic", ".seg")


## these will form the input files to gistic, after some post processing to get the headers in the right format
segments <- CNV.write(x, what = "segments")
segments_df <- bind_rows(segments)

data.table::fwrite(
  segments_df,
  file = segments_file,
  sep = "\t",
  col.names = TRUE
)

#write gistic input ?? 
gistic <- CNV.write(x, what="gistic", file=gistic_file)





