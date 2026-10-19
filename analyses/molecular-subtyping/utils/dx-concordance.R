#' Compare subtypes with the historical pathology diagnosis
#'
#' `dx_groups` maps each pathology_diagnosis to a group, and `family_compat`
#' lists, per subtype family, the groups that are acceptable (the first is the
#' primary group). Broad histology is deliberately not used, since it is a
#' manual lookup of pathology_diagnosis.

# a pathology_diagnosis can hold several values separated by ";"
dx_to_groups <- function(dx, dx_groups) {
  map_chr(dx, function(d) {
    if (is.na(d)) return(NA_character_)
    parts <- str_trim(str_split(d, ";")[[1]])
    grp <- dx_groups$dx_group[match(parts, dx_groups$pathology_diagnosis)]
    paste(sort(unique(replace_na(grp, "UNMAPPED"))), collapse = ";")
  })
}

# first matching family regex wins, top to bottom
subtype_compatible_groups <- function(subtype, family_compat) {
  map_chr(subtype, function(s) {
    hit <- family_compat %>% filter(str_detect(s, subtype_regex)) %>% slice(1)
    if (nrow(hit) == 0) NA_character_ else hit$compatible_dx_groups
  })
}

#' @return concordant, compatible, discordant, no_informative_dx or
#'   subtype_family_unmapped
dx_concordance <- function(subtype, dx_group, family_compat) {
  not_informative <- c("UNINFORMATIVE", "UNMAPPED")
  pmap_chr(
    list(subtype_compatible_groups(subtype, family_compat), dx_group),
    function(compat, group) {
      if (is.na(group) || group %in% not_informative) return("no_informative_dx")
      got <- setdiff(str_split(group, ";")[[1]], not_informative)
      if (length(got) == 0) return("no_informative_dx")
      if (is.na(compat)) return("subtype_family_unmapped")
      allowed <- str_split(compat, ";")[[1]]
      if (allowed[1] %in% got) "concordant"
      else if (any(got %in% allowed)) "compatible"
      else "discordant"
    }
  )
}

#' Tumor grade group from the free-text diagnosis
#'
#' Used when pathology_diagnosis is missing or uninformative ("Other").
#' @return "HGG" when the text names a high-grade term (anaplastic, WHO grade
#'   III/IV, high grade, glioblastoma, malignant), "LGG" for a low-grade term
#'   (WHO grade I/II, low grade), and NA when the text does not say.
free_text_grade_group <- function(text) {
  text <- str_to_lower(replace_na(text, ""))
  high <- "anaplastic|who grade (iii|3|iv|4)|grade (iii|iv|3|4)\\b|high.?grade|glioblastoma|gbm|malignant"
  low <- "who grade (i|ii|1|2)\\b|grade (i|ii|1|2)\\b|low.?grade"
  case_when(
    str_detect(text, high) ~ "HGG",
    str_detect(text, low) ~ "LGG",
    TRUE ~ NA_character_
  )
}
