suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(dplyr))

#this will convert all output files to qs2 AND delete .parquet files if run. 

# set up optparse options
option_list <- list(
  make_option(opt_str = "--output_dir", type = "character", default = NULL,
              help = "The directory containing the output .parquet files",
              metavar = "character"),
  make_option(opt_str = "--output_prefix", type = "character", default = NULL,
              help = "output file prefix",
              metavar = "character"),
  make_option(opt_str = "--manifest", type = "character", default = NULL,
              help = "file manifest",
              metavar = "character"),
  
  make_option(
    c("--FORCE"),
    action = "store_true",
    default = FALSE,
    help = "will force conversion to .qs2 even if manifest is greater than 6000 rows"
  )
  
)
opt <- parse_args(OptionParser(option_list = option_list))
out_dir <- opt$output_dir
out_pref <- opt$output_prefix
FORCE <- opt$FORCE

out_dir <- 'test-out'
out_pref <- 'test'
manifest <- 'controls_and_dicer_manifest.tsv'
FORCE <- FALSE

man_df <- read_tsv(manifest, show_col_types = FALSE)

n <- nrow(man_df)

if (isTRUE(FORCE) || n < 8000) {
  message("Continuing: rows = ", n,
          ", FORCE = ", isTRUE(FORCE), " Converting to .qs2")
  
} else {
  message(
    "Stopping: manifest has ", n, " rows (>= 8000) and FORCE is FALSE. Parquet format for large datasets preferred",
    "Use --FORCE to override."
  )
  quit(status = 1)
}


files <- list.files(out_dir, full.names = FALSE)

parquet_files <- files[
  grepl(paste0("^", out_pref), files) &
    grepl("\\.parquet$", files, ignore.case = TRUE)
]


for (f in parquet_files) {
  full_path <- file.path(out_dir, f)
  message("Processing: ", f)
  # Fully load the parquet file
  df <- open_dataset(full_path) %>%
    collect()
  # Build output filename (.qs2)
  out_file <- file.path(
    out_dir,
    sub("\\.parquet$", ".qs2", f)
  )
  # Save as qs2
  qs_save(df, out_file)
  # Clean up memory
  rm(df)
  gc()
}