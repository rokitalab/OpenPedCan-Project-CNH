################################################################################
# check-subtyping.R
#
# Checks of the molecular subtyping results and of the integrated histologies
# table, to run after run-molecular-subtyping.sh and run-subtyping-integrate.sh
# when new data come in. Each check passes or is flagged with a count and a few
# examples; the full table is written to results/subtyping_checks.tsv. The
# nuances behind the checks are in DECISIONS.md and README.md.
#
# Usage: Rscript --vanilla analyses/molecular-subtyping/utils/check-subtyping.R
################################################################################

suppressMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
results_dir <- file.path(module_dir, "results")
integrate_dir <- file.path(root_dir, "analyses", "molecular-subtyping-integrate")

checks <- list()
add_check <- function(area, check, flagged, n, examples = NA_character_) {
  checks[[length(checks) + 1]] <<- tibble(
    area = area, check = check, status = if_else(flagged, "FLAG", "pass"),
    n = n, examples = if_else(flagged, examples, NA_character_))
}
examples_of <- function(x, k = 6) {
  x <- unique(na.omit(as.character(x)))
  paste0(paste(head(x, k), collapse = "; "), if_else(length(x) > k, "; ...", ""))
}
read_result <- function(file, ...) {
  if (!file.exists(file)) return(NULL)
  read_tsv(file, show_col_types = FALSE, guess_max = 100000, ...)
}

m <- read_result(file.path(results_dir, "molecular_subtype.tsv"))
calls <- read_result(file.path(results_dir, "methylation_calls.tsv"))
hist <- read_result(file.path(integrate_dir, "results", "histologies.tsv"))
map <- read_result(file.path(integrate_dir, "results", "subtype_diagnosis_map.tsv"))
stopifnot(!is.null(m), !is.null(calls))

### Methylation inputs -------------------------------------------------------------

add_check("methylation", "calls come only from DKFZ v12.8 and NIH v3.1",
          !all(na.omit(calls$methyl_source) %in% c("dkfz_v12.8", "nih_v3.1")),
          sum(!na.omit(calls$methyl_source) %in% c("dkfz_v12.8", "nih_v3.1")),
          examples_of(setdiff(na.omit(calls$methyl_source), c("dkfz_v12.8", "nih_v3.1"))))

# old classifier versions must not be read anywhere in the module code
code_files <- list.files(module_dir, pattern = "\\.(R|Rmd)$", full.names = TRUE)
old_version <- map_dfr(code_files, function(f) {
  lines <- readLines(f, warn = FALSE)
  hit <- which(str_detect(lines, "dkfz_v12_methylation|dkfz_v11|NIH_v2|nih_v2") &
                 !str_detect(lines, "^\\s*#"))
  tibble(file = basename(f), line = hit)
})
add_check("methylation", "no code reads DKFZ v11, v12 or NIH v2 columns",
          nrow(old_version) > 0, nrow(old_version),
          examples_of(paste0(old_version$file, ":", old_version$line)))

gaps <- read_result(file.path(results_dir, "methylation_crosswalk_gaps.tsv"))
add_check("methylation", "every high-confidence class is in the crosswalk",
          !is.null(gaps) && nrow(gaps) > 0, if (is.null(gaps)) 0L else nrow(gaps),
          if (is.null(gaps)) NA_character_ else examples_of(paste0(gaps$call, " (", gaps$source, ", ", gaps$n_high_conf_biospecimens, ")")))
add_check("methylation", "no high-confidence call is unmapped",
          any(calls$unmapped_high_conf_call %in% TRUE),
          sum(calls$unmapped_high_conf_call %in% TRUE),
          examples_of(calls$methyl_class[calls$unmapped_high_conf_call %in% TRUE]))

if (!is.null(hist)) {
  bti_ids <- read_tsv(file.path(module_dir, "input", "bti_methylation_classifications.tsv"),
                      col_select = Bioassay_ID, show_col_types = FALSE)$Bioassay_ID
  # the TARGET arrays are not classified in the BTI table (not subtyped here)
  arrays <- hist %>% filter(experimental_strategy == "Methylation",
                            !cohort %in% "TARGET")
  missing_arrays <- setdiff(arrays$Kids_First_Biospecimen_ID, bti_ids)
  add_check("methylation", "every methylation array in the histologies is classified",
            length(missing_arrays) > 0, length(missing_arrays),
            examples_of(missing_arrays))
}

### Subtype calls ------------------------------------------------------------------

review <- m %>% filter(needs_review %in% TRUE | !is.na(subtype_review_flag))
add_check("subtype", "tumors flagged for review (diagnosis discordant, classifiers disagree, ...)",
          nrow(review) > 0, nrow(review),
          examples_of(paste0(review$sample_id, " (", review$subtype_review_flag, ")")))

discordant <- m %>% filter(dx_concordance %in% "discordant", subtype_source %in% "methylation",
                           is.na(diagnosis_review))
add_check("subtype", "methylation contradicts the diagnosis and is not yet reviewed (the subtype is kept; review the diagnosis; add reviewed ones to input/reviewed_diagnosis_discordance.tsv)",
          nrow(discordant) > 0, nrow(discordant),
          examples_of(paste0(discordant$sample_id, ": ", str_trunc(discordant$pathology_diagnosis, 24),
                             " -> ", discordant$molecular_subtype)))

non_tumor <- m %>% filter(str_detect(replace_na(molecular_subtype, ""),
                                     "^(Control|Non-tumor)|INFLAM"))
add_check("subtype", "control or non-tumor methylation classes on tumor samples",
          nrow(non_tumor) > 0, nrow(non_tumor),
          examples_of(paste0(non_tumor$sample_id, ": ", non_tumor$molecular_subtype)))

h3_outside <- m %>% filter(str_detect(replace_na(subtype_review_flag, ""), "H3 variant outside"))
add_check("subtype", "H3 K27/K28 variant outside DMG/DHG (a mix-up or a real H3 K27M tumor)",
          nrow(h3_outside) > 0, nrow(h3_outside), examples_of(h3_outside$sample_id))

cell_lines <- m %>% filter(str_detect(replace_na(subtype_review_flag, ""), "cell line differs"))
add_check("subtype", "cell line subtype differs from its tumor",
          nrow(cell_lines) > 0, nrow(cell_lines), examples_of(cell_lines$sample_id))

mycn <- m %>% filter(str_detect(replace_na(subtype_review_flag, ""), "MYCN"))
add_check("subtype", "neuroblastoma MYCN status differs between copy number, expression and the clinical report",
          nrow(mycn) > 0, nrow(mycn), examples_of(mycn$sample_id))

mismatch <- read_result(file.path(module_dir, "input", "possible_methylation_mixups.tsv"))
likely <- if (is.null(mismatch)) NULL else filter(mismatch, mismatch_call %in% "likely")
add_check("subtype", "likely methylation mismatches (genomics supports the diagnosis): collect and remove from the cohort (results/methylation_likely_mismatch.tsv)",
          !is.null(likely) && nrow(likely) > 0, if (is.null(likely)) 0L else nrow(likely),
          if (is.null(likely)) NA_character_ else examples_of(likely$sample_id))
possible <- if (is.null(mismatch)) NULL else filter(mismatch, mismatch_call %in% "possible")
add_check("subtype", "possible methylation mismatches still to investigate (input/possible_methylation_mixups.tsv)",
          !is.null(possible) && nrow(possible) > 0, if (is.null(possible)) 0L else nrow(possible),
          if (is.null(possible)) NA_character_ else examples_of(possible$sample_id))
# every tumor with a mismatch call is in the subtype table
if (!is.null(mismatch)) {
  # a likely mismatch has its array excluded (01-build-methylation-calls.R); a tumor with
  # no other evidence (no DNA, no free text rule) has no subtype row
  not_found <- setdiff(mismatch$sample_id[mismatch$mismatch_call != "likely"], m$sample_id)
  add_check("subtype", "every sample of the mismatch table is in the subtype table",
            length(not_found) > 0, length(not_found), examples_of(not_found))
}

# medulloblastoma: an SHH subtype from a low-scoring methylation call in an age
# group that does not fit it (gamma is an infant subtype)
shh <- read_result(file.path(results_dir, "mb_shh_subtypes_w_molecular_data.tsv"))
if (!is.null(shh)) {
  odd <- shh %>% filter(molecular_subtype %in% "MB, SHH gamma",
                        age_at_diagnosis_years >= 5,
                        shh_methylation_source %in% "dkfz_v12.8") %>%
    distinct(match_id, sample_id, age_at_diagnosis_years)
  add_check("subtype", "SHH gamma calls older than 5 years (check the methylation score and the UMAP)",
            nrow(odd) > 0, nrow(odd),
            examples_of(paste0(odd$sample_id, " (", round(odd$age_at_diagnosis_years), " y)")))
}

### Integrated table ---------------------------------------------------------------

if (!is.null(hist) && !is.null(map)) {
  tum <- hist %>% filter(sample_type == "Tumor")
  for (col in c("harmonized_diagnosis", "broad_histology", "short_histology")) {
    n <- sum(is.na(tum[[col]]))
    add_check("integrate", paste("every tumor has", col), n > 0, n,
              examples_of(tum$Kids_First_Biospecimen_ID[is.na(tum[[col]])]))
  }

  sub <- tum %>% filter(!is.na(molecular_subtype))
  assigned <- sub %>%
    filter(!str_detect(molecular_subtype, "To be classified|(?<!ETMR), NOS$|(?<!H3 )wildtype$"))
  n <- sum(is.na(assigned$integrated_diagnosis))
  add_check("integrate", "every subtyped tumor (specific subtype) has an integrated diagnosis",
            n > 0, n, examples_of(assigned$molecular_subtype[is.na(assigned$integrated_diagnosis)]))

  # "X, To be classified" (no molecular evidence) keeps the label; its integrated
  # diagnosis is the harmonized diagnosis
  placeholder <- sub %>% filter(str_detect(molecular_subtype, "To be classified"),
                                is.na(integrated_diagnosis))
  add_check("integrate", "'To be classified' labels have an integrated diagnosis (the harmonized diagnosis)",
            nrow(placeholder) > 0, n_distinct(placeholder$match_id),
            examples_of(placeholder$molecular_subtype))

  gaps <- map %>% filter(is.na(broad_histology) | is.na(short_histology))
  add_check("integrate", "subtype labels without a derived broad or short histology (the old diagnosis' value is kept; see results/subtype_histology_gaps.tsv)",
            nrow(gaps) > 0, nrow(gaps), examples_of(gaps$molecular_subtype))

  unmapped_cg <- read_result(file.path(integrate_dir, "results", "cancer_group_unmapped.tsv"))
  add_check("integrate", "every tumor diagnosis is in input/cancer_group_table.tsv (results/cancer_group_unmapped.tsv)",
            !is.null(unmapped_cg) && nrow(unmapped_cg) > 0,
            if (is.null(unmapped_cg)) 0L else nrow(unmapped_cg),
            if (is.null(unmapped_cg)) NA_character_ else examples_of(unmapped_cg$harmonized_diagnosis))

  plot_mapping_file <- file.path(root_dir, "figures", "palettes", "plot-mapping.tsv")
  if (file.exists(plot_mapping_file)) {
    plot_mapping <- read_tsv(plot_mapping_file, show_col_types = FALSE, na = character(),
                             col_types = cols(.default = col_character()))
    # cohorts that are not in the plot mapping (excluded_cohorts of extend-plot-mapping.R)
    not_plotted <- c("DGD", "TCGA", "PPTC", "Maris", "CPTAC")
    pairs <- tum %>% filter(!cohort %in% not_plotted) %>%
      mutate(cancer_group = replace_na(cancer_group, "NA")) %>%
      distinct(cohort, broad_histology, cancer_group) %>%
      anti_join(plot_mapping, by = c("cohort", "broad_histology", "cancer_group"))
    add_check("integrate", "every cohort, broad histology and cancer group has a plot group (figures/palettes/plot-mapping.tsv; add rows with figures/scripts/extend-plot-mapping.R)",
              nrow(pairs) > 0, nrow(pairs),
              examples_of(paste(pairs$cohort, pairs$broad_histology, pairs$cancer_group, sep = " / ")))
    to_review <- sum(plot_mapping$status != "")
    add_check("integrate", "plot-mapping rows added automatically and still to review (status column)",
              to_review > 0, to_review, NA_character_)
  }

  # a cancer group left empty is a deliberate row of input/cancer_group_table.tsv
  # (benign or non-tumor diagnoses); a diagnosis not in the table is checked above
  add_check("integrate", "every subtype label has a diagnosis rule",
            any(map$rule == "none"), sum(map$rule == "none"),
            examples_of(map$molecular_subtype[map$rule == "none"]))

  inconsistent <- sub %>% filter(!is.na(integrated_diagnosis)) %>%
    group_by(match_id) %>% filter(n_distinct(integrated_diagnosis) > 1) %>% ungroup()
  add_check("integrate", "all biospecimens of a tumor share the integrated diagnosis",
            nrow(inconsistent) > 0, n_distinct(inconsistent$match_id),
            examples_of(inconsistent$match_id))
}

### Reproducibility: the reviews -----------------------------------------------------
# The manual decisions are rows of tables keyed on sample_id (pathology reports and
# reviewer decisions are not in the repository), so a new release can leave a row
# without its tumor, make it redundant, or change the data it was made on.

read_input <- function(name) {
  f <- file.path(module_dir, "input", name)
  if (!file.exists(f)) return(NULL)
  read_tsv(f, show_col_types = FALSE, col_types = cols(.default = "c"))
}
review_tables <- list(
  manual_subtype_review.tsv = read_input("manual_subtype_review.tsv"),
  possible_methylation_mixups.tsv = read_input("possible_methylation_mixups.tsv"),
  reviewed_diagnosis_discordance.tsv = read_input("reviewed_diagnosis_discordance.tsv"),
  legacy_comparison_review.tsv = read_input("legacy_comparison_review.tsv"))
base_samples <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                         col_select = sample_id, show_col_types = FALSE)$sample_id

# every reviewed sample_id is a tumor of the current data (an orphan review is dead)
orphans <- imap_dfr(review_tables, function(x, name) {
  if (is.null(x)) return(NULL)
  tibble(table = name, sample_id = setdiff(x$sample_id, base_samples))
})
add_check("reviews", "every reviewed sample_id exists in the histologies",
          nrow(orphans) > 0, nrow(orphans),
          examples_of(paste0(orphans$table, ": ", orphans$sample_id)))

# every review has the person or report behind it and a note that says what was found
min_chars <- list(
  manual_subtype_review.tsv = c(source = 3, note = 15),
  possible_methylation_mixups.tsv = c(genomics_support = 15),
  reviewed_diagnosis_discordance.tsv = c(decision = 3, source = 3, note = 15),
  legacy_comparison_review.tsv = c(reviewed = 2, note = 15))
incomplete <- imap_dfr(review_tables, function(x, name) {
  if (is.null(x)) return(NULL)
  need <- min_chars[[name]]
  bad <- rep(FALSE, nrow(x))
  for (cc in names(need)) {
    bad <- bad | is.na(x[[cc]]) | nchar(str_trim(x[[cc]])) < need[[cc]]
  }
  tibble(table = name, sample_id = x$sample_id[bad])
})
add_check("reviews", "every review has a source and a note that says what was found",
          nrow(incomplete) > 0, nrow(incomplete),
          examples_of(paste0(incomplete$table, ": ", incomplete$sample_id)))

# a manual subtype is a label of the vocabulary of the module (a label that no tumor
# or rule uses is a typo)
manual <- review_tables$manual_subtype_review.tsv
known_labels <- unique(c(na.omit(m$molecular_subtype),
                         read_input("subtype_lesion_rules.tsv")$OPC_molecular_subtype))
unknown_label <- if (is.null(manual)) character() else
  setdiff(manual$molecular_subtype, known_labels)
add_check("reviews", "every manual subtype is a label used by the module",
          length(unknown_label) > 0, length(unknown_label), examples_of(unknown_label))

# stale: the automated subtype now equals the manual one (the row can be retired),
# or the methylation call a mix-up row was made on has changed
if (!is.null(manual) && "automated_subtype" %in% names(m)) {
  redundant <- m %>% filter(!is.na(manual_source), manual_conflict %in% FALSE) %>%
    distinct(sample_id)
  add_check("reviews", "manual reviews that the automated subtype now agrees with (can be retired)",
            nrow(redundant) > 0, nrow(redundant), examples_of(redundant$sample_id))
}
mixups <- review_tables$possible_methylation_mixups.tsv
if (!is.null(mixups)) {
  current_class <- calls %>% filter(!is.na(methyl_class)) %>%
    inner_join(read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                        col_select = c(Kids_First_Biospecimen_ID, sample_id),
                        show_col_types = FALSE), by = "Kids_First_Biospecimen_ID") %>%
    group_by(sample_id) %>%
    summarise(current = paste(unique(str_replace_all(methyl_class, "_", " ")), collapse = "; "),
              .groups = "drop")
  changed_class <- mixups %>%
    mutate(recorded = str_remove(str_replace_all(methyl_class, "_", " "), "\\s*\\(.*$")) %>%
    inner_join(current_class, by = "sample_id") %>%
    filter(!is.na(recorded), recorded != "", !str_detect(current, fixed(recorded))) %>%
    pull(sample_id)
  add_check("reviews", "mix-up rows whose methylation class differs from the current call",
            length(changed_class) > 0, length(changed_class), examples_of(changed_class))
}

# the data the committed results were made on (results/run_provenance.tsv, written
# by 07-compile-subtypes.R): a different input means the results are to be rerun
prov_file <- file.path(results_dir, "run_provenance.tsv")
if (file.exists(prov_file)) {
  prov <- read_tsv(prov_file, show_col_types = FALSE, col_types = cols(.default = "c"))
  now_md5 <- unname(tools::md5sum(file.path(root_dir, prov$file)))
  changed <- prov$file[!is.na(now_md5) & now_md5 != prov$md5]
  missing <- prov$file[is.na(now_md5)]
  add_check("reviews", "the input files are those of the committed results (results/run_provenance.tsv)",
            length(c(changed, missing)) > 0, length(c(changed, missing)),
            examples_of(c(changed, missing)))
} else {
  add_check("reviews", "the input files are those of the committed results (results/run_provenance.tsv)",
            TRUE, 1L, "results/run_provenance.tsv is missing: run 07-compile-subtypes.R")
}

### Output -------------------------------------------------------------------------

out <- bind_rows(checks)
write_tsv(out, file.path(results_dir, "subtyping_checks.tsv"))
print(out %>% mutate(check = str_trunc(check, 78), examples = str_trunc(examples, 60)),
      n = Inf, width = Inf)
message("\n", sum(out$status == "FLAG"), " of ", nrow(out),
        " checks flagged; see results/subtyping_checks.tsv")
