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
suppressWarnings(
  suppressPackageStartupMessages(library(tibble))
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
              'Bioassay_ID' columns. An optional 'sample_type' column excludes
              the samples it calls 'Normal' from tumor purity estimation."),
  make_option(
        opt_str = "--output_basename",
        type = "character", default = NULL,
        help = "The absolute path of the base directory containing sample array IDAT files.",
        metavar = "character"
    ),
  make_option(opt_str = "--funnorm", type = "logical", default = TRUE,
              help = "preprocesses the Illumina methylation arrays using one of
              the following minfi normalization methods: 
              - preprocessFunnorm: when array dataset contains either control
                                   samples (i.e., normal and tumor samples) or 
                                   multiple OpenPedcan cancer groups (TRUE)
              - preprocessQuantile: when an array dataset has only tumor samples
                                    from a single OpenPedcan cancer group (FALSE)
              Default is TRUE (preprocessFunnorm)",
              metavar = "TRUE|FALSE"),
  
  make_option(opt_str = "--snp_filter", type = "logical", default = TRUE,
              help = "If TRUE, drops the probes that contain either a SNP at
              the CpG interrogation or at the single nucleotide extension.
              Default is TRUE",
              metavar = "TRUE|FALSE"),
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


# read manifest to obtain the IDAT prefix from the `file_name` and its matched
# `Bioassay_ID` column, along with `sample_type` where the manifest carries it
man_df <- read_tsv(file = manifest_file, show_col_types = FALSE) %>% 
  dplyr::select(file_name, Bioassay_ID, dplyr::any_of("sample_type")) %>%
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
# force = TRUE is required: IDAT batches of the same array type can still
# report different array "sizes" to minfi (e.g. control-probe count differs
# across manufacturing lots), which read.metharray.exp() otherwise treats as
# an error. Confirmed 2026 (jharenza): removing force=TRUE fails on
# real EPIC batches with "different array size but seemingly all of the
# same type" even though 00-unzip-and-sort.R has already sorted by array type.
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

# Compute control-probe MAD separately for each channel. A zero-MAD channel
# indicates a degenerate readout even if the other channel has variation.
green_control_mad <- vapply(seq_len(ncol(green)), function(i) {
  mad(green[control_idx, i], na.rm = TRUE)
}, numeric(1))
red_control_mad <- vapply(seq_len(ncol(red)), function(i) {
  mad(red[control_idx, i], na.rm = TRUE)
}, numeric(1))
sample_names <- colnames(green)

# Free memory from large intermediate objects
rm(green, red, controls_info)
gc()

# Identify samples with a degenerate control-probe readout in either channel.
bad_sample_idx <- which(green_control_mad == 0 | red_control_mad == 0)
bad_samples <- sample_names[bad_sample_idx]

if (length(bad_samples) > 0) {
  message("Samples with MAD = 0 in at least one control channel (will be skipped):")
  zero_mad_samples <- paste0(out_base, "-", dataset, "-zero-mad.txt")
  fileConn <- file(zero_mad_samples)
  writeLines(bad_samples, fileConn)
  close(fileConn)
  print(bad_samples)

  # Apply this array-level QC filter before either normalization method.
  RGset <- RGset[, -bad_sample_idx, drop = FALSE]
  if (ncol(RGset) == 0) {
    stop("All samples have MAD = 0 in at least one control channel.")
  }

} else {
  message("No samples with MAD = 0 in either control channel detected.")
}

# Save the QC-filtered RGChannelSet for downstream CNV calling.
rg_set_file <- paste0(out_base, "-", dataset, "-rg-set.qs2")
qs_save(RGset, rg_set_file)

######################### Extract identity SNP probes #########################
# The arrays carry a small panel of rs probes that genotype common SNPs rather
# than measure methylation, for checking sample identity against genomic data.
# They are read straight off the RGChannelSet because they have no annotation
# entry, and so never reach the GenomicRatioSet the measurements below come
# from. This is unrelated to --snp_filter, which drops CpG probes overlapping a
# SNP. Betas are saved rather than genotype calls so the calling thresholds
# stay with the analysis making the comparison.
message("\nExtracting identity SNP probes...\n")

snp_file <- paste0(out_base, "-", dataset, "-methyl-snps.parquet")

snp_betas <- minfi::getSnpBeta(RGset) %>%
  as_tibble(rownames = "ProbeID") %>%
  rename_with(~ recode(.x, !!!setNames(man_df$Bioassay_ID, man_df$file_name)))

message("Identity SNP probes on ", dataset, ": ", nrow(snp_betas))

write_parquet(snp_betas, snp_file, compression = "zstd")

# minfi::detectionP() iterates over the Type I and Type II probe tables only,
# so it returns nothing at all for the rs probes. The background model it uses
# applies to them unchanged, because the Snp probes share the structure of the
# ordinary probes: SnpI probes measure both alleles in one channel, SnpII
# probes measure one address across both. This reproduces that calculation for
# the Snp addresses, and agrees with detectionP() to floating-point precision
# when run against the Type I and Type II probes instead.
# It mirrors the matrix method of minfi:::.detectionP(), which is where
# detectionP() does its arithmetic, with the Snp probe tables in place of the
# Type I and Type II ones. The upper tail is taken with lower.tail = FALSE
# rather than as 1 - pnorm(), which underflows to exactly 0 for a strong
# probe instead of returning the small p-value.
snp_detection_p <- function(RGset) {
  controls <- minfi::getProbeInfo(RGset, type = "Control")
  red <- minfi::getRed(RGset)
  green <- minfi::getGreen(RGset)
  negative <- intersect(controls$Address[controls$Type == "NEGATIVE"],
                        rownames(red))
  red_mu <- matrixStats::colMedians(red[negative, , drop = FALSE])
  red_sd <- matrixStats::colMads(red[negative, , drop = FALSE])
  green_mu <- matrixStats::colMedians(green[negative, , drop = FALSE])
  green_sd <- matrixStats::colMads(green[negative, , drop = FALSE])

  snp_i <- minfi::getProbeInfo(RGset, type = "SnpI")
  snp_ii <- minfi::getProbeInfo(RGset, type = "SnpII")
  i_red <- snp_i[snp_i$Color == "Red", ]
  i_green <- snp_i[snp_i$Color == "Grn", ]

  probe_names <- c(i_red$Name, i_green$Name, snp_ii$Name)
  detP <- matrix(NA_real_, nrow = length(probe_names), ncol = ncol(RGset),
                 dimnames = list(probe_names, colnames(RGset)))
  for (i in seq_len(ncol(RGset))) {
    if (nrow(i_red) > 0L) {
      intensity <- red[i_red$AddressA, i] + red[i_red$AddressB, i]
      detP[i_red$Name, i] <- pnorm(q = intensity, mean = red_mu[i] * 2,
                                   sd = red_sd[i] * 2, lower.tail = FALSE)
    }
    if (nrow(i_green) > 0L) {
      intensity <- green[i_green$AddressA, i] + green[i_green$AddressB, i]
      detP[i_green$Name, i] <- pnorm(q = intensity, mean = green_mu[i] * 2,
                                     sd = green_sd[i] * 2, lower.tail = FALSE)
    }
    if (nrow(snp_ii) > 0L) {
      intensity <- red[snp_ii$AddressA, i] + green[snp_ii$AddressA, i]
      detP[snp_ii$Name, i] <- pnorm(q = intensity,
                                    mean = red_mu[i] + green_mu[i],
                                    sd = red_sd[i] + green_sd[i],
                                    lower.tail = FALSE)
    }
  }
  detP
}

snp_p_file <- paste0(out_base, "-", dataset, "-methyl-snp-p-values.parquet")

snp_detP <- snp_detection_p(RGset) %>%
  as_tibble(rownames = "ProbeID") %>%
  rename_with(~ recode(.x, !!!setNames(man_df$Bioassay_ID, man_df$file_name)))

message("Identity SNP probe measurements with detection p > 0.01: ",
        sum(snp_detP[, -1] > 0.01, na.rm = TRUE), " of ",
        nrow(snp_detP) * (ncol(snp_detP) - 1))

write_parquet(snp_detP, snp_p_file, compression = "zstd")

rm(snp_betas, snp_detP)
gc()

########################## Estimate tumor purity ##############################
# RFpurify's ABSOLUTE and ESTIMATE random forests were trained on TCGA 450k
# betas computed with minfi's Illumina offset, so purity is predicted from a
# raw MethylSet rather than from the normalized values generated below.
# predict.randomForest() aborts on any NA in its input, so neither detection
# p-value masking nor SNP filtering is applied to the purity matrix.
#
# Both forests are defined on 450k probes, which 450k and EPICv1 carry in
# full. EPICv2 carries 1483 of the 1674 and is included anyway: blanking the
# missing 191 in 450k and EPICv1 samples, where the true value is known, moves
# ABSOLUTE by at most 0.014 and ESTIMATE by at most 0.043, at correlations of
# 0.998 and 0.983. n_probes_imputed records how much was filled per run.
purity_array_types <- c("IlluminaHumanMethylation450k",
                        "IlluminaHumanMethylationEPIC",
                        "IlluminaHumanMethylationEPICv2")
array_type <- minfi::annotation(RGset)[["array"]]

# Both forests were trained on tumors, to estimate what fraction of a tumor
# sample is tumor. A normal sample has no such fraction, but the forests still
# return a plausible-looking number for one rather than refusing, so normals
# are held out here and never reach the purity output. They are recognized by
# the manifest's optional `sample_type` column; without it there is nothing to
# exclude on and every sample is estimated, as before.
normal_samples <- character(0)
if ("sample_type" %in% names(man_df)) {
  normal_samples <- man_df$file_name[
    !is.na(man_df$sample_type) & tolower(man_df$sample_type) == "normal"
  ]
}
purity_samples <- setdiff(colnames(RGset), normal_samples)

if (!array_type %in% purity_array_types) {
  message("\nSkipping tumor purity estimation: RFpurify does not support ",
          array_type, ".\n")
} else if (!requireNamespace("RFpurify", quietly = TRUE)) {
  warning("RFpurify is not installed; skipping tumor purity estimation.")
} else if (length(purity_samples) == 0L) {
  message("\nSkipping tumor purity estimation: the manifest calls every ",
          dataset, " sample normal.\n")
} else {
  message("\nEstimating tumor purity with RFpurify...\n")

  if (!"sample_type" %in% names(man_df)) {
    warning("Manifest has no 'sample_type' column; estimating purity for ",
            "every sample, normal samples included.")
  }
  message("Samples for purity estimation: ", length(purity_samples), " of ",
          ncol(RGset), " (", ncol(RGset) - length(purity_samples),
          " normal excluded)")

  # The forests have to be attached, not merely loaded: predict_purity_betas()
  # refers to them by bare name, which does not reach the namespace's lazy-data
  # environment that `RFpurify::` resolves against. Attaching also masks
  # minfi::combine with randomForest::combine, which this script never calls.
  suppressWarnings(
    suppressPackageStartupMessages(library(RFpurify))
  )

  purity_file <- paste0(out_base, "-", dataset, "-methyl-purity.parquet")

  # The two forests use 856 probes each and share only 38 of them, so taking
  # the union lets one matrix serve both.
  model_probes <- unique(c(
    rownames(RFpurify::RFpurify_ABSOLUTE$importance),
    rownames(RFpurify::RFpurify_ESTIMATE$importance)
  ))

  # preprocessRaw() builds full Meth and Unmeth matrices covering every probe
  # on the array, of which only the model probes are wanted. It works one
  # sample at a time, so chunking over samples returns identical values while
  # keeping those transient matrices proportional to the chunk rather than to
  # the whole batch.
  purity_idx <- match(purity_samples, colnames(RGset))
  sample_chunks <- split(
    purity_idx,
    ceiling(seq_along(purity_idx) / 50)
  )
  # EPICv2 names every probe for its replicate, as cg00000029_TC21, so the
  # model probe names match nothing until the suffix is stripped. Stripping is
  # a no-op on 450k and EPICv1, whose locus names carry no suffix. 33 model
  # probes have more than one EPICv2 replicate; their betas are averaged,
  # since the replicates measure the same CpG and the forests take one value
  # per locus. The merge step picks a single replicate on detection p-value
  # instead, which is not available this early in the script.
  model_betas <- do.call(cbind, lapply(sample_chunks, function(idx) {
    m_set <- minfi::preprocessRaw(RGset[, idx])
    locus <- sub("_.*$", "", rownames(m_set))
    keep <- locus %in% model_probes
    betas <- minfi::getBeta(m_set[keep, ], type = "Illumina")
    rowsum(betas, locus[keep]) / as.vector(table(locus[keep]))
  }))
  gc()

  # predict_purity_betas() looks probes up by name and returns an all-NA row
  # for any it cannot find, which would abort the prediction. Give every model
  # probe a row, fill the ones this array lacks with 0.5, and report the count
  # so the estimates can be judged against how much was imputed.
  purity_betas <- matrix(
    NA_real_,
    nrow = length(model_probes), ncol = ncol(model_betas),
    dimnames = list(model_probes, colnames(model_betas))
  )
  purity_betas[rownames(model_betas), ] <- model_betas

  n_probes_imputed <- sum(!complete.cases(purity_betas))
  message("RFpurify model probes missing from ", array_type, ": ",
          n_probes_imputed, " of ", length(model_probes))
  purity_betas[is.na(purity_betas)] <- 0.5

  purity_values <- tibble(
    file_name = colnames(purity_betas),
    purity_absolute = as.numeric(
      RFpurify::predict_purity_betas(purity_betas, method = "ABSOLUTE")
    ),
    purity_estimate = as.numeric(
      RFpurify::predict_purity_betas(purity_betas, method = "ESTIMATE")
    ),
    array_type = array_type,
    n_model_probes = length(model_probes),
    n_probes_imputed = n_probes_imputed
  ) %>%
    mutate(Bioassay_ID = dplyr::recode(
      file_name, !!!setNames(man_df$Bioassay_ID, man_df$file_name)
    )) %>%
    select(Bioassay_ID, purity_absolute, purity_estimate, array_type,
           n_model_probes, n_probes_imputed)

  write_parquet(purity_values, purity_file, compression = "zstd")

  rm(model_betas, purity_betas, purity_values)
  gc()
}

######################## Calculate detection p-values #########################
message("\nsetting parallel processing options...\n")
library(BiocParallel)
BiocParallel::register(BiocParallel::MulticoreParam(workers = n_cores)) #UP THIS ON BIGGER MACHINE!!!
message("\nCalculating detection p-values...\n")




if (n_cores > 1) {
  message(sprintf("Running detectionP in parallel with %d cores...", n_cores))
  
  # split sample indices into chunks
  n <- ncol(RGset)
  chunks <- split(seq_len(n), ceiling(seq_len(n)/500))
  
  # Run detectionP in parallel on each chunk. Iterate over chunk positions
  # rather than over the chunks themselves, so the progress message reports the
  # chunk number instead of being vectorised over every sample index in it.
  det_list <- bplapply(seq_along(chunks), function(i) {
    message(sprintf("Processing chunk %d / %d", i, length(chunks)))
    minfi::detectionP(RGset[, chunks[[i]]])
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
sample_qc_file <- paste0(out_base, "-", dataset, "-methyl-sample-qc.parquet")


message("Extracting m values")

# extract m values (compute once, use for both masked and unmasked)
m_values <- minfi::getM(GRset)

m_values_unmasked <- m_values %>%
  as_tibble(rownames = "ProbeID") %>%
  rename_with(~ recode(.x, !!!setNames(man_df$Bioassay_ID, man_df$file_name)))

write_parquet(m_values_unmasked, m_value_file, compression = "zstd")

# Free memory
rm(m_values_unmasked)
gc()

##masking is optional for m values -- can generate masked and unmasked matrices

# Create masked version from the same m_values matrix
m_values[detP > 0.01] <- NA

m_values_masked <- m_values %>%
  as_tibble(rownames = "ProbeID") %>%
  rename_with(~ recode(.x, !!!setNames(man_df$Bioassay_ID, man_df$file_name)))

write_parquet(m_values_masked, m_value_file_masked, compression = "zstd")

# Free memory
rm(m_values, m_values_masked)
gc()

message("Extracting beta-values")

# Extract beta values and apply masking
beta_values <- minfi::getBeta(GRset)
beta_values[detP > 0.01] <- NA

beta_values_masked <- beta_values %>%
  as_tibble(rownames = "ProbeID") %>%
  rename_with(~ recode(.x, !!!setNames(man_df$Bioassay_ID, man_df$file_name)))
# Free beta_values matrix
rm(beta_values)
gc()
write_parquet(beta_values_masked, beta_value_file, compression = "zstd")

# ensure tibble
detP <- as_tibble(detP, rownames = "ProbeID")
# rename columns
colnames(detP) <- dplyr::recode(
  colnames(detP),
  !!!setNames(man_df$Bioassay_ID, man_df$file_name)
)

write_parquet(detP, p_value_file, compression = "zstd")

# Per-sample QC metric: fraction of probes failing detection (p > 0.01),
# for a downstream analysis-level sample-exclusion QC step.
sample_qc <- detP %>%
  select(-ProbeID) %>%
  summarise(across(everything(), ~ mean(.x > 0.01, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Bioassay_ID", values_to = "frac_failed_probes")

write_parquet(sample_qc, sample_qc_file, compression = "zstd")

# Free memory
rm(detP, beta_values_masked, sample_qc)
gc()

message("Extracting copy number values")
cn_value <- GRset %>% minfi::getCN() 
# ensure tibble
cn_value <- as_tibble(cn_value, rownames = "ProbeID")
# rename columns
colnames(cn_value) <- dplyr::recode(
  colnames(cn_value),
  !!!setNames(man_df$Bioassay_ID, man_df$file_name)
)


# write output file

write_parquet(cn_value, cn_value_file, compression = "zstd")
# delete GenomicRatioSet object to free memory
rm(GRset)
gc()
