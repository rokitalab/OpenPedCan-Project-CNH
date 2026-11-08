#' SQL source for a release table
#'
#' Returns the FROM source of a duckdb query: the parquet copy of the table in
#' data/parquet (from scripts/download-parquet.sh) when it exists, otherwise the
#' release TSV read as text. The columns have the same names in both.
#'
#' @param name Parquet table name (file name without .parquet).
#' @param tsv_file Release file name in data/ to fall back on.
#' @param root_dir Repository root.
release_source <- function(name, tsv_file, root_dir) {
  parquet <- file.path(root_dir, "data", "parquet", paste0(name, ".parquet"))
  if (file.exists(parquet)) {
    message("Reading ", name, " from parquet")
    sprintf("read_parquet('%s')", normalizePath(parquet))
  } else {
    message("Reading ", tsv_file, " (no parquet copy; run scripts/download-parquet.sh)")
    sprintf(paste0("read_csv('%s', delim = '\t', header = true, all_varchar = true,",
                   " sample_size = 100, ignore_errors = true)"),
            file.path(root_dir, "data", tsv_file))
  }
}
