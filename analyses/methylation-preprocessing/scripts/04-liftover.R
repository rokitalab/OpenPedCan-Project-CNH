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
gistic_seg_files <- gistic_seg_files[
  basename(gistic_seg_files) != "combined_hg38.gistic.seg"
]

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
seg_hg38 <- all_segs %>%
  filter(platform == "EPICv2") %>%
  mutate(Chromosome = as.character(Chromosome))
seg_hg19 <- all_segs %>% filter(platform %in% c("450k", "EPIC"))


resolve_liftover_overlaps <- function(df) {
  # Chain gaps can move adjacent hg19 segment boundaries past one another in
  # hg38.  GISTIC rejects overlapping segments for a sample/chromosome, so
  # split each overlap at its midpoint after ordering the lifted intervals.
  df <- df %>% arrange(Sample, Chromosome, Start_Position, End_Position)
  groups <- split(seq_len(nrow(df)), interaction(df$Sample, df$Chromosome, drop = TRUE))

  for (idx in groups) {
    if (length(idx) < 2) next

    for (i in 2:length(idx)) {
      previous <- idx[i - 1]
      current <- idx[i]

      if (df$Start_Position[current] <= df$End_Position[previous]) {
        overlap_start <- max(df$Start_Position[previous], df$Start_Position[current])
        overlap_end <- min(df$End_Position[previous], df$End_Position[current])
        boundary <- floor((overlap_start + overlap_end) / 2)
        df$End_Position[previous] <- boundary
        df$Start_Position[current] <- boundary + 1
      }
    }
  }

  df %>% filter(Start_Position <= End_Position)
}


lift_seg <- function(df, chain) {
  gr <- GRanges(
    seqnames = paste0("chr", df$Chromosome),
    ranges = IRanges(start = df$Start_Position, end = df$End_Position)
  )
  
  # Lift the two segment boundaries independently.  Lifting a full, often
  # chromosome-scale segment returns one range per chain block, which is not a
  # usable SEG representation for GISTIC.  A segment is retained when both
  # boundaries map uniquely to the same target chromosome.
  start_gr <- resize(gr, width = 1, fix = "start")
  end_gr <- resize(gr, width = 1, fix = "end")
  lifted_start <- liftOver(start_gr, chain)
  lifted_end <- liftOver(end_gr, chain)

  keep <- lengths(lifted_start) == 1 & lengths(lifted_end) == 1

  if (!any(keep)) {
    return(df[FALSE, ])
  }

  lifted_start <- unlist(lifted_start[keep])
  lifted_end <- unlist(lifted_end[keep])
  same_chromosome <- as.character(seqnames(lifted_start)) ==
    as.character(seqnames(lifted_end))

  if (!any(same_chromosome)) {
    return(df[FALSE, ])
  }

  df <- df[which(keep)[same_chromosome], , drop = FALSE]
  lifted_start <- lifted_start[same_chromosome]
  lifted_end <- lifted_end[same_chromosome]
  
  # Replace coordinates
  df$Chromosome <- as.character(seqnames(lifted_start))
  df$Chromosome <- gsub("^chr", "", df$Chromosome)
  df$Start_Position <- pmin(start(lifted_start), start(lifted_end))
  df$End_Position <- pmax(start(lifted_start), start(lifted_end))
  
  resolve_liftover_overlaps(df)
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
