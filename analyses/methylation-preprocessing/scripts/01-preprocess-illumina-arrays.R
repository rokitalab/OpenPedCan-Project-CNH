# Prepocess raw Illumina Infinium HumanMethylation BeadArrays (450K, and 850k) 
# intensities using minfi into usable methylation measurements (Beta and M values) 
# and copy number (cn-values) for OpenPedCan.

# Eric Wafula for Pediatric OpenTargets
# 09/28/2022

# Load libraries:
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(arrow))
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
  make_option(
        opt_str = "--output_basename",
        type = "character", default = NULL,
        help = "The absolute path of the base directory containing sample array IDAT files.",
        metavar = "character"
    ),
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
out_base <- opt$output_basename


#base_dir <- 'sorted_idats_output_dir/IlluminaHumanMethylationEPIC'
#snp_filter <- TRUE
#use_funnorm <- TRUE
#manifest_file <- 'controls_and_dicer_manifest.tsv'
#n_cores <- 4 
#out_base <- 'test-out/'


# read manifest to obtain the IDAT prefix from the `file_name` and its matched `Bioassay_ID` column
man_df <- read_tsv(file = manifest_file, show_col_types = FALSE) %>% 
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

# Compute MAD per sample by combining channels on the fly (avoids creating large stacked matrix)
control_mad <- sapply(seq_len(ncol(green)), function(i) {
  mad(c(green[control_idx, i], red[control_idx, i]), na.rm = TRUE)
})
names(control_mad) <- colnames(green)

# Free memory from large intermediate objects
rm(green, red, controls_info)
gc()

# Identify problematic samples
bad_samples <- names(control_mad[control_mad == 0])

if (use_funnorm) {
  
  if (length(bad_samples) > 0) {
    message("Samples with MAD = 0 (will be skipped):")
    zero_mad_samples <- paste0(out_base, "-", dataset, "zero-mad.txt")
    fileConn <- file(zero_mad_samples)
    writeLines(bad_samples, fileConn)
    close(fileConn)
    print(bad_samples)
    
    # Filter out bad samples
    RGset <- RGset[, control_mad > 0]
    
  } else {
    message("No samples with MAD = 0 detected.")
  }
  
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

# delete RGChannelSet object immediately to free memory
rm(RGset)
gc()

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
m_value_file <- paste0(out_base, "-", dataset, "-methyl-m-values-unmasked.parquet")
m_value_file_masked <- paste0(out_base, "-", dataset, "-methyl-m-values-masked.parquet")
beta_value_file <- paste0(out_base, "-", dataset, "-methyl-beta-values-masked.parquet")
cn_value_file <- paste0(out_base, "-", dataset, "-methyl-cn-values.parquet")
p_value_file <- paste0(out_base, "-", dataset, "-methyl-p-values.parquet")
m_set_file <- paste0(out_base, "-", dataset, '-m-set.qs2')


qs_save(GRset, m_set_file)

message("Extracting m values")

# extract m values (compute once, use for both masked and unmasked)
m_values <- minfi::getM(GRset)

# Create unmasked version
m_values_unmasked <- m_values %>% as.data.frame() %>%
  tibble::rownames_to_column("Probe_ID")

m_values_unmasked <- data.table::setnames(m_values_unmasked, man_df$file_name, man_df$Bioassay_ID, skip_absent = TRUE)

# write output file

if (!is.data.frame(m_values_unmasked)) {
  m_values_masked <- as.data.frame(m_values_unmasked)
}
write_parquet(m_values_unmasked, m_value_file)
#qs_save(m_value_unmasked, m_value_file)

# Free memory
rm(m_values_unmasked)
gc()

##masking is optional for m values -- can generate masked and unmasked matrices

# Create masked version from the same m_values matrix
m_values[detP > 0.05] <- NA
m_value_masked <- m_values %>% as.data.frame() %>%
  tibble::rownames_to_column("Probe_ID")

m_value_masked <- data.table::setnames(m_value_masked, man_df$file_name, man_df$Bioassay_ID, skip_absent = TRUE)

# write output file
if (!is.data.frame(m_value_masked)) {
  m_value_masked <- as.data.frame(m_value_masked)
}
write_parquet(m_value_masked, m_value_file_masked)

#qs_save(m_value_masked, m_value_file_masked)

# Free memory
rm(m_values, m_value_masked)
gc()

message("Extracting beta-values")

# Extract beta values and apply masking
beta_values <- minfi::getBeta(GRset)
beta_values[detP > 0.05] <- NA
beta_values_masked <- beta_values %>% as.data.frame() %>%
  tibble::rownames_to_column("Probe_ID")

# Free beta_values matrix
rm(beta_values)
gc()

beta_values_masked <- data.table::setnames(beta_values_masked, man_df$file_name, man_df$Bioassay_ID, skip_absent = TRUE)

# write output file
if (!is.data.frame(beta_values_masked)) {
  beta_values_masked <- as.data.frame(beta_values_masked)
}
write_parquet(beta_values_masked, beta_value_file)
#qs_save(beta_values_masked, beta_value_file)

if (!is.data.frame(detP)) {
  detP <- as.data.frame(detP)
}
detP <- data.table::setnames(detP, man_df$file_name, man_df$Bioassay_ID, skip_absent = TRUE)
# write output file

write_parquet(detP, p_value_file)

# Free memory
rm(detP, beta_values_masked)
gc()

message("Extracting copy number values")
cn_value <- GRset %>% minfi::getCN() %>% as.data.frame() #%>% #keep as tibble if saving as qs2 ? 
  #tibble::rownames_to_column("Probe_ID")
if (!is.data.frame(cn_value)) {
  cn_value <- as.data.frame(cn_value)
}
cn_value <- data.table::setnames(cn_value, man_df$file_name, man_df$Bioassay_ID, skip_absent = TRUE)

# write output file

write_parquet(cn_value, cn_value_file)
#qs_save(cn_value, cn_value_file)
# delete GenomicRatioSet object to free memory
rm(GRset)
gc()
