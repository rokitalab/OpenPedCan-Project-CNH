#' Open an in-memory duckdb connection
#'
#' duckdb >= 1.5.6 prints a note about its temporary extension directory unless
#' `shared_home` is set; older versions do not have that argument, so it is only
#' passed when it exists.
connect_duckdb <- function() {
  args <- list()
  if ("shared_home" %in% names(formals(duckdb::duckdb))) {
    args$shared_home <- FALSE
  }
  DBI::dbConnect(do.call(duckdb::duckdb, args))
}
