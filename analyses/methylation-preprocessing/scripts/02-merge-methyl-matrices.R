# Merges methylation beta-values, m-values, and cp-values matrices for 
# all pre-processed array datasets.

# Jessica Daggett
# 05/12/2026

# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(minfi))

# Magrittr pipe
`%>%` <- dplyr::`%>%`


# set up optparse options
option_list <- list(
  make_option(opt_str = "--output_dir", type = "character", default = NULL,
              help = "The directory containing the output .parquet files",
              metavar = "character")
)

# establish base dir
#root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

out_dir <- 'test-out'

#find out what array types are present
files <- list.files(out_dir, full.names = FALSE)
array_types <- character()

if (any(grepl("EPICv2", files, ignore.case = TRUE))) {
  array_types <- c(array_types, "EPICv2")
}
if (
  any(grepl("EPICv1", files, ignore.case = TRUE)) ||
  any(grepl("EPIC", files, ignore.case = TRUE) &
      !grepl("EPICv2", files, ignore.case = TRUE))
) {
  array_types <- c(array_types, "EPICv1")
}
if (any(grepl("450k", files, ignore.case = TRUE))) {
  array_types <- c(array_types, "450k")
}

detP_epicv2_fn <- file.path(out_dir, "EPICv2-IlluminaHumanMethylationEPICv2-methyl-p-values.parquet")

detP <- open_dataset(detP_epicv2_fn)

detP_probe <- detP %>% ##summarises detP to get median p val for each probe
  group_by(probeID) %>%
  summarise(
    detP_median = median(detP, na.rm = TRUE),
    .groups = "drop"
  )

v2_probes <- epicv2 %>%
  select(probeID) %>%
  collect() %>%
  mutate(Probe_base = sub("_.*$", "", probeID)) ## collects probes and gets 'probe base' with _BCxxx _TCxxx removed

v2_probes_unique <- v2_probes %>% ##uses the lowest detP median
  left_join(detP_probe, by = "probeID") %>%
  filter(Probe_base %in% v1_probes$probeID) %>%
  filter(!is.na(detP_median)) %>%
  group_by(Probe_base) %>%
  arrange(detP_median, probeID, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()
