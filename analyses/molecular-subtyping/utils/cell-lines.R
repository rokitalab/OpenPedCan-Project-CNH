#' Compare derived cell lines with their tumors
#'
#' A derived cell line is subtyped from its own methylation and lesions (its own
#' match_id) and compared with the solid-tissue tumors of the same sample (or,
#' failing that, of the same participant): `tumor_subtype`, `tumor_group`,
#' `tumor_fidelity` (same, same group, different, no subtype to compare) and
#' `tumor_alterations_shared` (the tumor's confirming alterations that the line
#' carries, lacks or could not show). A line with no subtype of its own that
#' shares a confirming lesion with its tumor takes the tumor subtype, and a line
#' that still differs is flagged `cell line differs from tumor`. Also writes
#' results/cell_line_tumor_fidelity.tsv.
#'
#' @param compiled The compiled subtype table (one row per match_id).
#' @param root_dir Repository root.
#' @param results_dir The module results directory (it holds
#'   cell_line_alterations.tsv from 03-validate-subtypes.R).
#' @return `compiled` with the cell line columns added.
add_cell_line_comparison <- function(compiled, root_dir, results_dir) {
  ### Cell lines versus the tumor --------------------------------------
  # A derived cell line is subtyped from its own methylation and lesions (its own
  # match_id). The tumors are the solid-tissue tumors of the same sample
  # (or, failing that, of the same participant). Concordance: same = identical
  # subtype, same group = the same umbrella group (e.g. HGG, H3 wildtype for an
  # HGG, RTK tumor), different = another subtype, and no subtype on either side.
  composition <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                          guess_max = 100000, show_col_types = FALSE) %>%
    filter(!is.na(match_id)) %>%
    distinct(match_id, composition, cohort)
  # every derived line is kept out of the tumors; only the PBTA lines are
  # compared with a parental tumor (the neuroblastoma lines have no primary tumor)
  is_line_all <- composition$match_id[composition$composition %in% "Derived Cell Line"]
  is_line <- composition$match_id[composition$composition %in% "Derived Cell Line" &
                                    composition$cohort %in% "PBTA"]

  tumors_by_line <- compiled %>%
    filter(!match_id %in% is_line_all, !is.na(molecular_subtype),
           !is.na(sample_id), !is.na(Kids_First_Participant_ID)) %>%
    select(Kids_First_Participant_ID, sample_id, tumor_subtype = molecular_subtype,
           tumor_group = subtype_group)
  by_sample <- tumors_by_line %>% group_by(sample_id) %>%
    summarise(tumor_subtype = paste(sort(unique(tumor_subtype)), collapse = ";"),
              tumor_group = paste(sort(unique(tumor_group)), collapse = ";"),
              .groups = "drop")
  by_participant <- tumors_by_line %>% group_by(Kids_First_Participant_ID) %>%
    summarise(tumor_subtype_p = paste(sort(unique(tumor_subtype)), collapse = ";"),
              tumor_group_p = paste(sort(unique(tumor_group)), collapse = ";"),
              .groups = "drop")

  compiled <- compiled %>%
    mutate(is_cell_line = match_id %in% is_line) %>%
    left_join(by_sample, by = "sample_id", na_matches = "never") %>%
    left_join(by_participant, by = "Kids_First_Participant_ID", na_matches = "never") %>%
    mutate(
      tumor_subtype = if_else(is_cell_line, coalesce(tumor_subtype, tumor_subtype_p),
                               NA_character_),
      tumor_group = if_else(is_cell_line, coalesce(tumor_group, tumor_group_p),
                             NA_character_),
      tumor_fidelity = case_when(
        !is_cell_line ~ NA_character_,
        is.na(molecular_subtype) | is.na(tumor_subtype) ~ "no subtype to compare",
        molecular_subtype == tumor_subtype ~ "same",
        subtype_group == tumor_group ~ "same group",
        TRUE ~ "different"
      )
    ) %>%
    select(-tumor_subtype_p, -tumor_group_p, -is_cell_line)

  ### Alterations the line shares with its tumor ----------------------------------
  # The alterations that confirm the tumor's subtype, checked in the line:
  # shared, not found (the assay was available) or not assessable (no DNA / RNA).
  # A line with no subtype of its own (no methylation, no lesion rule) that shares
  # a lesion with its tumor takes the tumor's subtype.
  line_alterations <- read_tsv(file.path(results_dir, "cell_line_alterations.tsv"),
                               show_col_types = FALSE)

  # the alterations of the tumors, with the subtype they confirm
  tumor_alts <- compiled %>%
    filter(!match_id %in% is_line_all, !is.na(subtype_alterations_confirmed),
           !is.na(sample_id)) %>%
    select(sample_id, Kids_First_Participant_ID, tumor_subtype = molecular_subtype,
           alts = subtype_alterations_confirmed) %>%
    separate_rows(alts, sep = "; ") %>%
    filter(!str_detect(alts, "^expression classifier"), alts != "") %>%
    distinct()

  assess_line <- function(sample_id, participant, line_subtype, line_source, items,
                          has_dna, has_rna) {
    wanted <- tumor_alts %>%
      filter(sample_id == !!sample_id |
               (is.na(!!sample_id) & Kids_First_Participant_ID == !!participant))
    if (nrow(wanted) == 0) wanted <- tumor_alts %>%
      filter(Kids_First_Participant_ID == !!participant)
    if (nrow(wanted) == 0) return(list(text = NA_character_, adopt = NA_character_,
                                       shared = NA_character_))
    have <- if (is.na(items)) character() else str_split(items, "; ")[[1]]
    rna_item <- function(x) str_detect(x, "--|overexpression")
    status <- map_chr(wanted$alts, function(a) {
      case_when(
        a %in% have ~ "shared",
        rna_item(a) & !has_rna ~ "not assessable (no RNA)",
        !rna_item(a) & !has_dna ~ "not assessable (no DNA)",
        TRUE ~ "not found in line"
      )
    })
    w <- wanted %>% mutate(status = status)
    parts <- w %>% group_by(status) %>%
      summarise(items = paste(sort(unique(alts)), collapse = ", "), .groups = "drop") %>%
      mutate(text = paste0(status, ": ", items))
    shared <- w %>% filter(status == "shared") %>% pull(alts) %>% unique()
    adopt <- NA_character_
    subtypes <- unique(w$tumor_subtype[w$status == "shared"])
    # not the umbrella wildtype or an annotation-only item such as TP53
    if (length(subtypes) == 1 && length(shared) > 0 &&
        (is.na(line_subtype) || line_source %in% "lesion negative")) {
      adopt <- subtypes
      # EZHIP is part of the label only when the line overexpresses it
      if (str_detect(adopt, ", EZHIP") && !"EZHIP overexpression" %in% shared) {
        adopt <- str_remove(adopt, ", EZHIP")
      }
    }
    list(text = paste(parts$text, collapse = "; "), adopt = adopt,
         shared = paste(sort(shared), collapse = "; "))
  }

  line_rows <- compiled %>%
    filter(!is.na(tumor_fidelity)) %>%
    left_join(line_alterations, by = "match_id") %>%
    mutate(res = pmap(
      list(sample_id, Kids_First_Participant_ID, molecular_subtype,
           subtype_source, line_alterations, has_dna, has_rna),
      assess_line
    )) %>%
    transmute(match_id,
              tumor_alterations_shared = map_chr(res, "text"),
              adopt_subtype = map_chr(res, "adopt"),
              adopt_alterations = map_chr(res, "shared"))

  compiled <- compiled %>%
    left_join(line_rows, by = "match_id") %>%
    mutate(
      adopted = !is.na(adopt_subtype),
      molecular_subtype = if_else(adopted, adopt_subtype, molecular_subtype),
      subtype_group = if_else(adopted, adopt_subtype, subtype_group),
      subtype_source = if_else(adopted, "lesion shared with tumor", subtype_source),
      evidence_level = if_else(adopted, "lesion shared with tumor", evidence_level),
      subtype_alterations_confirmed = if_else(adopted, adopt_alterations,
                                              subtype_alterations_confirmed),
      alteration_confirmed = alteration_confirmed | adopted,
      Notes = if_else(adopted, "Subtype from lesion shared with the tumor", Notes),
      tumor_fidelity = case_when(
        is.na(tumor_fidelity) ~ NA_character_,
        !is.na(molecular_subtype) & molecular_subtype == tumor_subtype ~ "same",
        # one label is the other plus an annotation term (e.g. ", EGFR", "-mutant")
        !is.na(molecular_subtype) & !is.na(tumor_subtype) &
          (str_starts(tumor_subtype, fixed(molecular_subtype)) |
             str_starts(molecular_subtype, fixed(tumor_subtype))) ~ "same group",
        adopted & subtype_group == tumor_group ~ "same group",
        TRUE ~ tumor_fidelity
      )
    ) %>%
    select(-adopt_subtype, -adopt_alterations, -adopted)
  message("Cell lines taking the tumor subtype from a shared lesion: ",
          sum(compiled$subtype_source %in% "lesion shared with tumor"))

  # a line that still differs from its tumor is flagged for review
  compiled <- compiled %>%
    mutate(
      needs_review = needs_review | tumor_fidelity %in% "different",
      subtype_review_flag = case_when(
        tumor_fidelity %in% "different" & is.na(subtype_review_flag) ~
          "cell line differs from tumor",
        tumor_fidelity %in% "different" ~
          paste(subtype_review_flag, "cell line differs from tumor", sep = ";"),
        TRUE ~ subtype_review_flag
      )
    )

  write_tsv(
    compiled %>% filter(!is.na(tumor_fidelity)) %>%
      select(match_id, sample_id, molecular_subtype, subtype_source,
             tumor_subtype, tumor_fidelity, tumor_alterations_shared),
    file.path(results_dir, "cell_line_tumor_fidelity.tsv")
  )

  compiled
}
