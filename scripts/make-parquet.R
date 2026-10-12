################################################################################
# make-parquet.R
#
# Write parquet copies of the large OpenPedCan release tables, for querying
# with duckdb / duckplyr / arrow (locally or from S3) without reading the whole
# TSV or RDS file. The TSV / RDS files stay the source of truth; the parquet
# tables are generated from them and checked against them (row counts).
#
# Output: data/<release>/parquet/<table>.parquet, plus
#   manifest.tsv         table, source file, rows, columns, size, md5
#   md5sum-parquet.txt   md5 and file name, in the md5sum.txt format
#
# Tables are sorted by gene (or position) with small row groups, so a gene
# lookup reads only the matching row groups. Column types are detected over the
# whole file (sample_size = -1), not from the first rows.
#
# Usage: Rscript scripts/make-parquet.R [release]   (default v15)
#   Set PARQUET_TABLES to a comma-separated list of table names to run a subset.
################################################################################

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

args <- commandArgs(trailingOnly = TRUE)
release <- if (length(args) >= 1) args[1] else "v15"

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data", release)
out_dir <- file.path(data_dir, "parquet")
tmp_dir <- file.path(root_dir, "scratch", "parquet-duckdb-tmp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

# table name, source file, ORDER BY (SQL), row group size, and optionally the
# body of a SELECT * REPLACE (...) for a documented value normalization
tsv_tables <- list(
  list("snv-consensus-plus-hotspots", "snv-consensus-plus-hotspots.maf.tsv.gz",
       "Hugo_Symbol, Chromosome, Start_Position", 100000),
  list("snv-mutect2-tumor-only-plus-hotspots",
       "snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz",
       "Hugo_Symbol, Chromosome, Start_Position", 100000),
  # v15 mixes Gain / Loss / Neutral (Control-FREEC tumor-only input) with
  # lowercase calls in status; lowercased here, and fixed upstream in
  # analyses/focal-cn-file-preparation/04-prepare-cn-file.R
  list("consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only",
       "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz",
       "gene_symbol, biospecimen_id", 100000, "lower(status) AS status"),
  list("consensus_seg_with_status", "consensus_seg_with_status.tsv",
       "Kids_First_Biospecimen_ID, chrom, \"loc.start\"", 100000),
  list("cnvkit_with_status", "cnvkit_with_status.tsv",
       "Kids_First_Biospecimen_ID, chrom, \"loc.start\"", 100000),
  list("fusion-annoFuse", "fusion-annoFuse.tsv.gz", "FusionName, Sample", 50000),
  list("fusion-arriba", "fusion-arriba.tsv.gz", "tumor_id, gene1, gene2", 50000),
  list("fusion-starfusion", "fusion-starfusion.tsv.gz", "tumor_id, FusionName", 50000),
  list("fusion-dgd", "fusion-dgd.tsv.gz", "FusionName, Sample", 50000),
  list("fusion-putative-oncogenic", "fusion-putative-oncogenic.tsv",
       "FusionName, Sample", 50000),
  list("sv-manta", "sv-manta.tsv.gz", "SV_chrom, SV_start", 100000)
)

# genes x biospecimens matrices: a gene_symbol column and one column per
# biospecimen, sorted by gene, with small row groups
rds_tables <- list(
  list("gene-expression-rsem-tpm-collapsed", "gene-expression-rsem-tpm-collapsed.rds"),
  list("gene-counts-rsem-expected_count-collapsed",
       "gene-counts-rsem-expected_count-collapsed.rds")
)

only <- Sys.getenv("PARQUET_TABLES", "")
keep <- function(name) only == "" || name %in% strsplit(only, ",")[[1]]

# duckdb >= 1.5.6 prints a note about its temporary extension directory unless
# shared_home is set; older versions do not have the argument
duckdb_args <- list(dbdir = ":memory:")
if ("shared_home" %in% names(formals(duckdb::duckdb))) {
  duckdb_args$shared_home <- FALSE
}
con <- dbConnect(do.call(duckdb::duckdb, duckdb_args))
invisible(dbExecute(con, "SET memory_limit = '10GB'"))
invisible(dbExecute(con, "SET threads = 4"))
invisible(dbExecute(con, sprintf("SET temp_directory = '%s'", tmp_dir)))

manifest <- list()

record <- function(name, source, parquet, source_rows = NA_integer_) {
  rows <- dbGetQuery(con, sprintf("SELECT count(*) n FROM read_parquet('%s')",
                                  parquet))$n
  cols <- nrow(dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM read_parquet('%s')",
                                       parquet)))
  manifest[[name]] <<- data.frame(
    table = name, source = source, rows = rows, source_rows = source_rows,
    columns = cols, size_mb = round(file.size(parquet) / 1e6, 1),
    md5 = unname(tools::md5sum(parquet)), stringsAsFactors = FALSE
  )
  message(sprintf("  %s: %s rows, %d columns, %.1f MB", name,
                  format(rows, big.mark = ","), cols, file.size(parquet) / 1e6))
}

for (t in tsv_tables) {
  name <- t[[1]]
  if (!keep(name)) next
  src <- file.path(data_dir, t[[2]])
  out <- file.path(out_dir, paste0(name, ".parquet"))
  if (!file.exists(src)) { message("Skipping ", name, ": ", src, " not found"); next }
  message("Converting ", name)
  replace_sql <- if (length(t) >= 5) sprintf(" REPLACE (%s)", t[[5]]) else ""
  invisible(dbExecute(con, sprintf(
    "COPY (SELECT *%s FROM read_csv('%s', delim = '\t', header = true, quote = '',
                                  sample_size = -1, nullstr = ['NA', ''],
                                  strict_mode = false)
           ORDER BY %s)
     TO '%s' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE %d)",
    replace_sql, src, t[[3]], out, t[[4]])))
  n_src <- dbGetQuery(con, sprintf(
    "SELECT count(*) n FROM read_csv('%s', delim = '\t', header = true, quote = '',
                                     strict_mode = false)", src))$n
  record(name, t[[2]], out, n_src)
}

for (t in rds_tables) {
  name <- t[[1]]
  if (!keep(name)) next
  src <- file.path(data_dir, t[[2]])
  out <- file.path(out_dir, paste0(name, ".parquet"))
  if (!file.exists(src)) { message("Skipping ", name, ": ", src, " not found"); next }
  message("Converting ", name)
  mat <- readRDS(src)
  df <- data.frame(gene_symbol = rownames(mat), mat, check.names = FALSE,
                   stringsAsFactors = FALSE)
  rownames(df) <- NULL
  df <- df[order(df$gene_symbol), ]
  duckdb_register(con, "mat", df)
  invisible(dbExecute(con, sprintf(
    "COPY (SELECT * FROM mat ORDER BY gene_symbol)
     TO '%s' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 2000)", out)))
  duckdb_unregister(con, "mat")
  record(name, t[[2]], out, nrow(mat))
  rm(mat, df); invisible(gc())
}

manifest_df <- do.call(rbind, manifest)
manifest_df$row_count_matches <- manifest_df$rows == manifest_df$source_rows
manifest_file <- file.path(out_dir, "manifest.tsv")
if (only != "" && file.exists(manifest_file)) {
  old <- read.delim(manifest_file, stringsAsFactors = FALSE)
  manifest_df <- rbind(old[!old$table %in% manifest_df$table, names(manifest_df)],
                       manifest_df)
}
write.table(manifest_df, file.path(out_dir, "manifest.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
files <- list.files(out_dir, pattern = "\\.parquet$", full.names = TRUE)
writeLines(paste(unname(tools::md5sum(files)), basename(files), sep = "  "),
           file.path(out_dir, "md5sum-parquet.txt"))
dbDisconnect(con, shutdown = TRUE)
message("\nDone. ", nrow(manifest_df), " tables in ", out_dir)
print(manifest_df[, c("table", "rows", "source_rows", "columns", "size_mb",
                      "row_count_matches")], row.names = FALSE)
