# Prepocess raw Illumina Infinium HumanMethylation BeadArrays (450K, and 850k) 
# intensities using minfi into usable methylation measurements (Beta and M values) 
# and copy number (cn-values) for OpenPedCan.

# Eric Wafula for Pediatric OpenTargets
# 09/28/2022

# Load libraries:
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(qs2))
suppressWarnings(
  suppressPackageStartupMessages(library(minfi))
)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# set up optparse options
option_list <- list(
  make_option(opt_str = "--base_dir", type = "character", default = NULL,
              help = "The absolute path of the base directory containing sample 
              array IDAT files.",
              metavar = "character"),
  make_option(opt_str = "--manifest_file", type = "character",
              help = "Input manifest file with 'file_name' and
              'Bioassay_ID' columns"),
  make_option(opt_str = "--funnorm", action = "store_true", 
              default = TRUE,
              help = "preprocesses the Illumina methylation arrays using one of
              the following minfi normalization methods: 
              - preprocessFunnorm: when array dataset contains either control
                                   samples (i.e., normal and tumor samples) or 
                                   multiple OpenPedcan cancer groups (TRUE)
              - preprocessQuantile: when an array dataset has only tumor samples
                                    from a single OpenPedcan cancer group (FALSE)
              Default is TRUE (preprocessFunnorm)",
              metavar = "character"),
  
  make_option(opt_str = "--snp_filter", action = "store_true", default = TRUE, 
              help = "If TRUE, drops the probes that contain either a SNP at
              the CpG interrogation or at the single nucleotide extension.
              Default is TRUE",
              metavar = "character"),
  make_option(opt_str = "--n_cores", type = 'integer',
              default=1, help="number of cores for parallelisation of minfi::detectionP. Default is 1")
)


# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
base_dir <- opt$base_dir
use_funnorm <- opt$funnorm
snp_filter <- opt$snp_filter
manifest_file <- opt$manifest_file
n_cores <- opt$n_cores


#base_dir <- 'inputs'
#snp_filter <- TRUE
#use_funnorm <- TRUE
#manifest_file <- 'inputs/epicv2-test.tsv'
#n_cores <- 4 

# read manifest to obtain the IDAT prefix from the `file_name` and its matched `Bioassay_ID` column
man_df <- read_tsv(file = manifest_file) %>% 
  select(file_name, Bioassay_ID) %>%
  dplyr::mutate(file_name = gsub("(_Red|_Grn).*", "", file_name)) %>%
  dplyr::mutate(file_name = basename(file_name)) %>%
  unique()

# get analysis cancer type from arrays base_dir
dataset <- basename(base_dir)
message("===============================================")
message(c("Preprocessing ", dataset, " sample array data files..."))
message("===============================================\n")

########################### Read sample array data  ############################
message("Reading sample array data files...\n")

# load array data into a RGChannelSet object
RGset <- suppressWarnings(
  minfi::read.metharray.exp(base = base_dir, verbose = TRUE, force = TRUE, recursive = TRUE)
)
###################### Check for MAD=0 samples and skip them ####################
message("\nChecking for samples with zero MAD in control probes...\n")





# Extract raw intensities from RGset
green <- minfi::getGreen(RGset)
red   <- minfi::getRed(RGset)

#get control probes
controls_info <- minfi::getProbeInfo(RGset, type = "Control")
control_idx <- rownames(green) %in% controls_info$Address
# Combine both channels into one matrix for MAD calculation

ctrl_matrix <- rbind(green[control_idx, ], red[control_idx, ])
# Compute MAD per sample (column)
control_mad <- apply(ctrl_matrix, 2, mad, na.rm = TRUE)

# Identify problematic samples
bad_samples <- names(control_mad[control_mad == 0])

if (use_funnorm) {
  
  if (length(bad_samples) > 0) {
    message("Samples with MAD = 0 (will be skipped):")
    print(bad_samples)
    
    # Optionally show file paths
    if (!is.null(minfi::getPaths(RGset))) {
      message("\nAssociated IDAT file paths:")
      print(basename(minfi::getPaths(RGset)[bad_samples]))
    }
    
    # Filter out bad samples
    RGset <- RGset[, control_mad > 0]
    
  } else {
    message("No samples with MAD = 0 detected.")
  }
  
}

####################### Pre-processing and normalization ########################
message("\nPre-processing and normalizing...\n")

# process data into a GenomicRatioSet object
if (use_funnorm) {
  # preprocessFunnorm
  GRset <- RGset %>% 
    minfi::preprocessFunnorm(nPCs=2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE,
                             keepCN = TRUE, ratioConvert = TRUE, verbose = TRUE)
} else { 
  # processQuantile
  GRset <- RGset %>%  
    minfi::preprocessQuantile(fixOutliers = TRUE,  quantileNormalize = TRUE,
                              stratified = TRUE, mergeManifest = TRUE, sex = NULL)
}

######################## Calculate detection p-values #########################
message("\nsetting parallel processing options...\n")
library(BiocParallel)
BiocParallel::register(BiocParallel::MulticoreParam(workers = n_cores)) #UP THIS ON BIGGER MACHINE!!!
message("\nCalculating detection p-values...\n")




if (n_cores > 1) {
  message(sprintf("Running detectionP in parallel with %d cores...", n_cores))
  
  # register parallel backend
  register(MulticoreParam(workers = n_cores))
  
  # split sample indices into chunks
  n <- ncol(RGset)
  chunks <- split(seq_len(n), ceiling(seq_len(n)/500))
  
  # run detectionP in parallel on each chunk
  det_list <- bplapply(chunks, function(idx) {
    message(sprintf("Processing chunk %d / %d ", idx, length(chunks)))
    minfi::detectionP(RGset[, idx])
  }, BPPARAM = MulticoreParam(n_cores))
  
  # combine results
  detP <- do.call(cbind, det_list)
  
  # reset to single-core for safety
  register(SerialParam())
  
} else {
  message("Running detectionP serially on 1 core...")
  detP <- minfi::detectionP(RGset)
}


register(SerialParam())

# delete RGChannelSet object to free memory
rm(RGset)

if (snp_filter) {
  ########################## Remove probes with SNPs ############################
  message("\nRemoving probes with SNPs...\n")
  
  # removing probes with SNPs inside the probe body
  # or at the nucleotide extension
  GRset <- GRset %>% 
    minfi::addSnpInfo() %>% 
    minfi::dropLociWithSnps(snps=c("SBE","CpG"), maf=0)
}


# Match detP to probes remaining in GRset -- add this after the GRset is filtered for snps? 
detP <- detP[featureNames(GRset), ]



############################## Generate results ###############################
message("Generate results...\n")

# extract relevant methylation values, copy number values and probe annotations
# from the GenomicRatioSet object

# set up output file names
m_value_file <- paste0(dataset, "-methyl-m-values-unmasked.qs2")
m_value_file_masked <- paste0(dataset, "-methyl-m-values-masked.qs2")
beta_value_file <- paste0(dataset, "-methyl-beta-values-masked.qs2")
cn_value_file <- paste0(dataset, "-methyl-cn-values.qs2")

message("Extracting m values")

# extract m values
m_value_unmasked <- GRset %>% minfi::getM() %>% as.data.frame() %>%
  tibble::rownames_to_column("Probe_ID")

m_value_unmasked <- data.table::setnames(m_value_unmasked, man_df$file_name, man_df$Bioassay_ID, skip_absent = TRUE)




# write output file

qs_save(m_value_unmasked, m_value_file)
##masking is optional for m values -- can generate masked and unmasked matrices

# extract m values
m_value_masked <- GRset %>% minfi::getM() %>% { .[detP > 0.05] <- NA; . } %>% as.data.frame() %>%
  tibble::rownames_to_column("Probe_ID")

m_value_masked <- data.table::setnames(m_value_masked, man_df$file_name, man_df$Bioassay_ID, skip_absent = TRUE)

# write output file

qs_save(m_value_masked, m_value_file_masked)
message("Extracting beta-values")

#beta_value <- GRset %>% minfi::getBeta() %>% as.data.frame() %>%
  #tibble::rownames_to_column("Probe_ID")


beta_values_masked <- GRset %>% 
  minfi::getBeta() %>%
  { .[detP > 0.05] <- NA; . } %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Probe_ID")

# apply masking -- #should ALWAYS be done for B values 
#beta_values_masked <- beta_values
#beta_values_masked[detP > 0.05] <- NA
beta_values_masked <- data.table::setnames(beta_values_masked, man_df$file_name, man_df$Bioassay_ID, skip_absent = TRUE)

# write output file

qs_save(beta_values_masked, beta_value_file)
message("Extracting copy number values")
cn_value <- GRset %>% minfi::getCN() %>% as.data.frame() %>%
  tibble::rownames_to_column("Probe_ID")

cn_value <- data.table::setnames(cn_value, man_df$file_name, man_df$Bioassay_ID, skip_absent = TRUE)

# write output file

qs_save(cn_value, cn_value_file)
# delete GenomicRatioSet object to free memory
rm(GRset)
