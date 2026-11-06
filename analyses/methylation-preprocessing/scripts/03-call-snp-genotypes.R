# Call genotypes from the identity SNP probe betas written by
# 01-preprocess-illumina-arrays.R, for checking sample identity against
# genomic data.
# Jessica Daggett, 11/05/2026

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(arrow))

option_list <- list(
  make_option("--output_dir", type = "character", default = NULL,
              help = "Directory containing input and output files"),
  make_option("--output_prefix", type = "character", default = NULL,
              help = "Output file prefix"),
  make_option("--input_prefix", type = "character", default = NULL,
              help = "Comma-separated input file prefixes; defaults to output_prefix"),
  make_option("--het_min", type = "double", default = 0.25,
              help = "Betas at or below this are called homozygous. Default 0.25"),
  make_option("--het_max", type = "double", default = 0.75,
              help = "Betas above this are called homozygous. Default 0.75"),
  make_option("--min_detp", type = "double", default = 0.01,
              help = "Measurements with a detection p-value above this are left
              uncalled. Set to 1 to call every measurement. Default 0.01"),
  make_option("--vcf", type = "logical", default = FALSE,
              help = "Also write a VCF per array type. Default FALSE")
)
opt <- parse_args(OptionParser(option_list = option_list))
out_dir <- opt$output_dir
out_pref <- opt$output_prefix
if (is.null(out_dir) || !dir.exists(out_dir)) {
  stop("--output_dir must name an existing directory.")
}
if (is.null(out_pref) || !nzchar(trimws(out_pref))) {
  stop("--output_prefix is required.")
}
if (!is.finite(opt$het_min) || !is.finite(opt$het_max) ||
    opt$het_min < 0 || opt$het_max > 1 || opt$het_min >= opt$het_max) {
  stop("--het_min and --het_max must satisfy 0 <= het_min < het_max <= 1.")
}
if (!is.finite(opt$min_detp) || opt$min_detp <= 0 || opt$min_detp > 1) {
  stop("--min_detp must satisfy 0 < min_detp <= 1.")
}
mask_detp <- opt$min_detp < 1
input_prefixes <- if (is.null(opt$input_prefix)) {
  out_pref
} else {
  unique(trimws(strsplit(opt$input_prefix, ",", fixed = TRUE)[[1]]))
}
if (length(input_prefixes) == 0L || any(!nzchar(input_prefixes))) {
  stop("--input_prefix must contain nonempty prefixes.")
}

array_names <- c(
  EPICv1 = "IlluminaHumanMethylationEPIC",
  EPICv2 = "IlluminaHumanMethylationEPICv2",
  `450k` = "IlluminaHumanMethylation450k"
)

# Match exact filenames, as the merge step does, so that prefix text cannot
# pull in an unintended array type or batch.
input_files <- bind_rows(lapply(input_prefixes, function(prefix) {
  tibble(
    array_type = names(array_names),
    beta_path = file.path(out_dir, paste0(
      prefix, "-", array_names, "-methyl-snps.parquet")),
    detp_path = file.path(out_dir, paste0(
      prefix, "-", array_names, "-methyl-snp-p-values.parquet"))
  )
})) %>%
  filter(file.exists(beta_path))
if (nrow(input_files) == 0L) {
  stop("No SNP probe files found for prefixes '",
       paste(input_prefixes, collapse = ", "), "' in ", out_dir, ".")
}
# Preprocessing runs predating the SNP detection p-value output have betas but
# no p-values. Rather than quietly calling every measurement, say so and let
# the caller ask for that explicitly with --min_detp 1.
if (mask_detp && any(!file.exists(input_files$detp_path))) {
  stop("Missing SNP detection p-value files: ",
       paste(basename(input_files$detp_path[!file.exists(input_files$detp_path)]),
             collapse = ", "),
       ". Rerun preprocessing, or pass --min_detp 1 to call without masking.")
}
message("SNP probe files: ", nrow(input_files))

# The three genotype classes sit at beta ~0, ~0.5 and ~1. Which homozygote is
# reference and which is alternate is not recoverable from the array, so the
# 0/0 and 1/1 labels are arbitrary with respect to a VCF; resolving that
# orientation per probe is the job of the analysis doing the comparison.
# Betas outside [0, 1] are treated as uncalled rather than clamped.
call_genotype <- function(beta, het_min, het_max) {
  gt <- rep("./.", length(beta))
  callable <- !is.na(beta) & beta >= 0 & beta <= 1
  gt[callable & beta <= het_min] <- "0/0"
  gt[callable & beta > het_min & beta <= het_max] <- "0/1"
  gt[callable & beta > het_max] <- "1/1"
  gt
}

read_long <- function(path, value_column) {
  values <- read_parquet(path)
  if (!"ProbeID" %in% names(values) || ncol(values) < 2L) {
    stop("Expected ProbeID and at least one sample column in ", basename(path))
  }
  if (anyDuplicated(values$ProbeID)) {
    stop("Duplicate ProbeID values in ", basename(path))
  }
  pivot_longer(values, -ProbeID, names_to = "sample", values_to = value_column)
}

read_snp_calls <- function(beta_path, detp_path, array_type) {
  long <- read_long(beta_path, "beta")
  if (mask_detp) {
    long <- long %>%
      inner_join(read_long(detp_path, "detp"), by = c("ProbeID", "sample"))
    if (nrow(long) == 0L) {
      stop("No probes shared between ", basename(beta_path), " and ",
           basename(detp_path))
    }
  } else {
    long$detp <- NA_real_
  }
  long %>%
    # EPICv2 carries several probes for some sites, distinguished by a suffix.
    # rs_id is what joins the panel across array types.
    mutate(probe_id = ProbeID,
           rs_id = sub("_.*$", "", ProbeID),
           array_type = array_type,
           beta = round(beta, 4),
           gt = call_genotype(beta, opt$het_min, opt$het_max),
           # A measurement that failed detection carries no usable genotype,
           # and a missing p-value cannot show that it passed.
           gt = if (mask_detp) {
             ifelse(is.na(detp) | detp > opt$min_detp, "./.", gt)
           } else {
             gt
           }) %>%
    select(probe_id, rs_id, sample, array_type, beta, detp, gt)
}

genotypes <- bind_rows(lapply(seq_len(nrow(input_files)), function(i) {
  message("Reading ", basename(input_files$beta_path[i]))
  read_snp_calls(input_files$beta_path[i], input_files$detp_path[i],
                 input_files$array_type[i])
}))

duplicated_samples <- genotypes %>%
  distinct(sample, array_type) %>%
  count(sample) %>%
  filter(n > 1L)
if (nrow(duplicated_samples) > 0L) {
  stop("Sample names occur on more than one array type: ",
       paste(duplicated_samples$sample, collapse = ", "))
}

genotype_file <- file.path(out_dir,
                           paste0(out_pref, "-methyl-snp-genotypes.tsv"))
write_tsv(genotypes, genotype_file)
message("Saved: ", genotype_file)

n_called <- sum(genotypes$gt != "./.")
message(sprintf("%d probe-sample pairs, %d called (%.1f%%)",
                nrow(genotypes), n_called,
                100 * n_called / nrow(genotypes)))
if (mask_detp) {
  message(sprintf("Left uncalled for detection p > %s: %d",
                  opt$min_detp,
                  sum(is.na(genotypes$detp) | genotypes$detp > opt$min_detp)))
}
genotypes %>%
  count(array_type, gt) %>%
  pivot_wider(names_from = gt, values_from = n, values_fill = 0L) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# A VCF is written per array type rather than one combined file, because the
# array types interrogate different probe sets and a combined file would have
# to invent missing calls for every probe absent from an array.
#
# CHROM and POS are placeholders: the rs probes carry no coordinates in the
# array annotation, and REF and ALT are unknown for the reason given above.
# These files are therefore only usable by matching on the ID column, and will
# not support position-based or allele-based operations.
if (opt$vcf) {
  for (array_type in unique(genotypes$array_type)) {
    rows <- genotypes %>%
      filter(array_type == !!array_type) %>%
      arrange(probe_id, sample)
    samples <- sort(unique(rows$sample))
    # DETP is emitted only when the p-values were read, so that a run with
    # --min_detp 1 does not carry a column of missing values. %.3g keeps the
    # exponent on the small p-values a passing probe produces.
    format_field <- if (mask_detp) "GT:BT:DETP" else "GT:BT"
    wide <- rows %>%
      mutate(field = paste0(gt, ":",
                            ifelse(is.na(beta), ".", sprintf("%.4f", beta)),
                            if (mask_detp) {
                              paste0(":", ifelse(is.na(detp), ".",
                                                 sprintf("%.3g", detp)))
                            } else {
                              ""
                            })) %>%
      select(probe_id, sample, field) %>%
      pivot_wider(names_from = sample, values_from = field) %>%
      arrange(probe_id)

    vcf_file <- file.path(out_dir, paste0(
      out_pref, "-", array_names[[array_type]], "-methyl-snp-genotypes.vcf"
    ))
    header <- c(
      "##fileformat=VCFv4.2",
      "##source=03-call-snp-genotypes.R",
      "##NOTE=CHROM/POS are placeholders and REF/ALT are unknown. Match on ID only.",
      sprintf("##NOTE=Genotypes called at het_min=%s, het_max=%s, min_detp=%s",
              opt$het_min, opt$het_max, opt$min_detp),
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      '##FORMAT=<ID=BT,Number=1,Type=Float,Description="Methylation beta">',
      if (mask_detp) {
        '##FORMAT=<ID=DETP,Number=1,Type=Float,Description="Detection p-value; GT is ./. above min_detp">'
      },
      paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
              "FORMAT", samples), collapse = "\t")
    )
    body <- paste(
      ".", seq_len(nrow(wide)), wide$probe_id, "N", "A", ".", ".", ".",
      format_field,
      apply(as.matrix(wide[, samples, drop = FALSE]), 1, paste, collapse = "\t"),
      sep = "\t"
    )
    writeLines(c(header, body), vcf_file)
    message("Saved: ", vcf_file)
  }
}
