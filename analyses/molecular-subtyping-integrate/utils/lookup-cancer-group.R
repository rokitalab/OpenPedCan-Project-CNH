#' Cancer group of a harmonized diagnosis, from a table
#'
#' `input/cancer_group_table.tsv` has one row per harmonized diagnosis
#' (`harmonized_diagnosis`, exact match, for the few diagnoses that differ from
#' their token) or per diagnosis token (`diagnosis_token`: the text before the
#' first comma, without the part in brackets). An exact row wins over a token row.
#' An empty `cancer_group` is deliberate: the diagnosis has no cancer group
#' (benign, non-tumor, ...). A diagnosis in neither is unmapped: its cancer group
#' is missing and it is listed for review.
#'
#' @param harmonized_diagnosis Character vector of harmonized diagnoses.
#' @param table The cancer group table (harmonized_diagnosis, diagnosis_token,
#'   cancer_group).
#' @return A data frame with `cancer_group` and `mapped` (TRUE when a row was
#'   found, also for an empty cancer group).
lookup_cancer_group <- function(harmonized_diagnosis, table) {
  token <- stringr::str_trim(stringr::str_extract(
    stringr::str_extract(harmonized_diagnosis, "[^,]*"), "[^(]*"))
  exact <- table[!is.na(table$harmonized_diagnosis), ]
  by_token <- table[!is.na(table$diagnosis_token), ]
  stopifnot(!anyDuplicated(exact$harmonized_diagnosis),
            !anyDuplicated(by_token$diagnosis_token))
  i_exact <- match(harmonized_diagnosis, exact$harmonized_diagnosis)
  i_token <- match(token, by_token$diagnosis_token)
  mapped <- !is.na(i_exact) | !is.na(i_token)
  group <- dplyr::if_else(!is.na(i_exact), exact$cancer_group[i_exact],
                          by_token$cancer_group[i_token])
  data.frame(cancer_group = dplyr::na_if(group, ""), mapped = mapped,
             stringsAsFactors = FALSE)
}
