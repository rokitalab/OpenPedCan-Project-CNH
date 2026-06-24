library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# set up optparse options
option_list <- list(
  make_option(opt_str = "--seg_dir", type = "character", default = NULL,
              help = "The path where the seg files are located",
              metavar = "character")
)


# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
seg_dir <- opt$seg_dir

# Magrittr pipe
`%>%` <- dplyr::`%>%`


# Locate SEG files
gistic_seg_files <- list.files(seg_dir, pattern = "\\gistic.seg$", full.names = TRUE)


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

# -----------------------------
# 3. Split by platform
# -----------------------------
seg_epicv2 <- all_segs %>% filter(platform == "EPICv2")
seg_legacy <- all_segs %>% filter(platform %in% c("450k", "EPIC"))


# Liftover (hg19 to hg38)


# Load chain file
chain <- import.chain("liftover/hg19ToHg38.over.chain")
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


if (nrow(seg_legacy) > 0) {
  seg_legacy_hg38 <- lift_seg(seg_legacy, chain)
} else {
  seg_legacy_hg38 <- NULL
}


# Combine all segments (hg38)

seg_epicv2 <- seg_epicv2 %>%
  mutate(Chromosome = as.character(Chromosome))

seg_legacy_hg38 <- seg_legacy_hg38 %>%
  mutate(Chromosome = as.character(Chromosome))
combined_seg <- bind_rows(
  seg_epicv2,     # already hg38
  seg_legacy_hg38 # lifted
)

combined_seg <- combined_seg %>%
  dplyr::select(-platform, -source_file)

#Write output
write_tsv(combined_seg, file = file.path(seg_dir, "combined_hg38.gistic.seg"))