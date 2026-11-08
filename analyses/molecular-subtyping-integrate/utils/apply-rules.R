#' Apply an ordered rule table to a data frame
#'
#' Each rule sets `target` to `result` for the rows that meet all of its
#' conditions, when no earlier rule of the same stage and target has set the
#' row (first match wins, like `case_when()`). Rules of a stage are evaluated in
#' `order`; stages in increasing order. The conditions read the current values,
#' so a later stage sees the result of the earlier ones.
#'
#' Condition columns (a blank cell is no condition):
#' * a column named like a data frame column (e.g. `broad_histology`): the
#'   value equals one of the `|`-separated alternatives; a leading `!` negates
#'   the whole cell;
#' * `<column>_regex`: `grepl()` of the pattern on that column;
#' * `<column>_na`: `TRUE` when the column is NA.
#' `result` is the new value, `NA` for a missing value, or `=<column>` to copy
#' another column.
#'
#' @param df A data frame.
#' @param rules The rule table: stage, order, target, the condition columns and
#'   result.
#' @param stages The stages to apply (default: all, in increasing order).
#' @return `df` with the targets updated.
apply_rules <- function(df, rules, stages = sort(unique(rules$stage))) {
  condition_cols <- setdiff(names(rules), c("stage", "order", "target", "result", "note"))
  for (st in sort(stages)) {
    stage_rules <- rules[rules$stage == st, ] %>% dplyr::arrange(order)
    done <- list()
    for (i in seq_len(nrow(stage_rules))) {
      r <- stage_rules[i, ]
      target <- r$target
      if (is.null(done[[target]])) done[[target]] <- rep(FALSE, nrow(df))
      hit <- !done[[target]]
      for (col in condition_cols) {
        cell <- r[[col]]
        if (is.na(cell) || cell == "") next
        if (stringr::str_ends(col, "_regex")) {
          column <- stringr::str_remove(col, "_regex$")
          hit <- hit & grepl(cell, df[[column]]) & !is.na(df[[column]])
        } else if (stringr::str_ends(col, "_na")) {
          column <- stringr::str_remove(col, "_na$")
          hit <- hit & (is.na(df[[column]]) == as.logical(cell))
        } else {
          negate <- stringr::str_starts(cell, "!")
          values <- stringr::str_split(stringr::str_remove(cell, "^!"), "\\|")[[1]]
          match <- df[[col]] %in% values
          hit <- hit & (if (negate) !match & !is.na(df[[col]]) else match)
        }
      }
      if (!any(hit)) next
      value <- r$result
      new <- if (is.na(value) || value == "NA") {
        rep(NA_character_, sum(hit))
      } else if (stringr::str_starts(value, "=")) {
        as.character(df[[stringr::str_remove(value, "^=")]][hit])
      } else {
        value
      }
      df[[target]] <- as.character(df[[target]])
      df[[target]][hit] <- new
      done[[target]] <- done[[target]] | hit
    }
  }
  df
}
