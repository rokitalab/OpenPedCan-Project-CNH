library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(optparse)

# set up optparse options
option_list <- list(
  make_option(opt_str = "--seg_dir", type = "character", default = NULL,
              help = "The path where the seg files are located",
              metavar = "character"),
  make_option(opt_str = "--chain_file", type = "character", default = NULL,
              help = "hg19-to-hg38 chain file. Required when 450k or EPICv1 SEG files are present.",
              metavar = "character")
)


# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
seg_dir <- opt$seg_dir
chain_file <- opt$chain_file

# Magrittr pipe
`%>%` <- dplyr::`%>%`


# Locate SEG files
gistic_seg_files <- list.files(seg_dir, pattern = "gistic\\.seg$", full.names = TRUE)

if (length(gistic_seg_files) == 0) {
  stop("No GISTIC SEG files found in: ", seg_dir)
}


#Detect platform from filename

detect_platform <- function(f) {
  f <- basename(f)
  if (grepl("EPICv2", f, ignore.case = TRUE)) return("EPICv2")
  if (grepl("EPIC-", f, ignore.case = TRUE)) return("EPIC")
  if (grepl("450", f, ignore.case = TRUE)) return("450k")
  return("unknown")
}

seg_df_list <- lapply(gistic_seg_files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  df$platform <- detect_platform(f)
  df$source_file <- basename(f)
  df
})

all_segs <- bind_rows(seg_df_list)

if (any(all_segs$platform == "unknown")) {
  stop(
    "Could not determine array platform for: ",
    paste(unique(all_segs$source_file[all_segs$platform == "unknown"]), collapse = ", ")
  )
}

# EPICv2 annotations use hg38; EPICv1 and 450k annotations use hg19.
seg_hg38 <- all_segs %>% filter(platform == "EPICv2")
seg_hg19 <- all_segs %>% filter(platform %in% c("450k", "EPIC"))


lift_seg <- function(df, chain) {
  gr <- GRanges(
    seqnames = paste0("chr", df$Chromosome),
    ranges = IRanges(start = df$Start_Position, end = df$End_Position)
  )
  
  lifted <- liftOver(gr, chain)
  
  # Keep only segments that map uniquely to ONE region
  keep <- lengths(lifted) == 1
  
  gr_lifted <- unlist(lifted[keep])
  df <- df[keep, ]
  
  # Replace coordinates
  df$Chromosome <- as.character(seqnames(gr_lifted))
  df$Chromosome <- gsub("^chr", "", df$Chromosome)
  df$Start_Position <- start(gr_lifted)
  df$End_Position <- end(gr_lifted)
  
  df
}

# Standardize all legacy arrays to hg38 before producing the combined SEG.
# This is required for mixed EPICv2 input and keeps a 450k/EPICv1-only run
# ready for the same hg38 GISTIC workflow.
if (nrow(seg_hg19) > 0) {
  if (is.null(chain_file) || !file.exists(chain_file)) {
    stop("A valid --chain_file is required when 450k or EPICv1 SEG files are present.")
  }

  message("Lifting 450k/EPICv1 segments from hg19 to hg38.")
  chain <- import.chain(chain_file)
  combined_seg <- bind_rows(seg_hg38, lift_seg(seg_hg19, chain))
  output_file <- file.path(seg_dir, "combined_hg38.gistic.seg")
} else {
  message("Only EPICv2 (hg38) SEG files detected; no liftover required.")
  combined_seg <- seg_hg38
  output_file <- file.path(seg_dir, "combined_hg38.gistic.seg")
}

combined_seg <- combined_seg %>%
  mutate(Chromosome = as.character(Chromosome)) %>%
  dplyr::select(-platform, -source_file)

#Write output
write_tsv(combined_seg, file = output_file)
