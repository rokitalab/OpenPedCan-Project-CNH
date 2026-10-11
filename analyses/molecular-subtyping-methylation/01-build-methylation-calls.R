################################################################################
# 01-build-methylation-calls.R
#
# Build one resolved, per-biospecimen methylation call table that downstream
# molecular subtyping modules consume instead of re-implementing thresholds and
# classifier-version aliasing.
#
# Calls come from two sources:
#   - DKFZ v12.8 / NIH v3.1: input/bti_methylation_classifications.tsv
#   - DKFZ v12 / NIH v2:     histologies-base.tsv
#
# Confidence rule (same as molecular-subtyping-Rare-CNS):
#   high  = DKFZ subclass score >= 0.8, or (DKFZ missing / < 0.8) NIH
#           superfamily and class scores >= 0.9
# Precedence among passing calls: DKFZ v12.8, NIH v3.1, DKFZ v12, NIH v2
#
# Class vocabularies differ between classifier versions, so every call is
# mapped to a single crosswalk row via normalized keys instead of hand aliases.
#
# Usage: Rscript --vanilla 01-build-methylation-calls.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-methylation")
input_dir <- file.path(module_dir, "input")
results_dir <- file.path(module_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

dkfz_cutoff <- 0.8
nih_cutoff <- 0.9

# normalize vocabulary differences between versions (case, "_" vs " ")
norm_key <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("[ _]+", " ", x)
  x[x %in% c("", "NA")] <- NA_character_
  x
}

### Crosswalk -> long key table -------------------------------------------------

crosswalk <- read_excel(
  file.path(input_dir, "mnp_v12.5_v12.8_NIHv2_v3.1_annotation_with_OPC_subtypes.xlsx")
) %>%
  mutate(crosswalk_row = row_number())

key_cols <- c("dkfz_Abbreviation", "dkfz_Abbreviation_v12.8",
              "Abbreviation_internal_NIHv2", "Abbreviation_internal_NIHv3",
              "OPC_molecular_subtype_methyl")

key_table <- crosswalk %>%
  select(crosswalk_row, all_of(key_cols)) %>%
  pivot_longer(-crosswalk_row, values_to = "raw_key") %>%
  mutate(key = norm_key(raw_key)) %>%
  filter(!is.na(key)) %>%
  distinct(key, crosswalk_row)

# a key must resolve to exactly one crosswalk row, otherwise calls are ambiguous
ambiguous_keys <- key_table %>%
  count(key) %>%
  filter(n > 1)
if (nrow(ambiguous_keys) > 0) {
  warning("Ambiguous crosswalk keys dropped: ",
          paste(ambiguous_keys$key, collapse = ", "))
  key_table <- key_table %>% filter(!key %in% ambiguous_keys$key)
}

# the crosswalk marks deliberate no-call classes (control tissues, inflammatory
# microenvironment) with the literal string "NA" in OPC_molecular_subtype; a
# blank cell instead means the class has not been assigned a subtype yet
crosswalk_info <- crosswalk %>%
  transmute(
    crosswalk_row,
    no_call_class = OPC_molecular_subtype %in% "NA",
    OPC_molecular_subtype = if_else(no_call_class, NA_character_,
                                    OPC_molecular_subtype),
    rare = Rare %in% "Yes",
    molecular_superfamily = Molecular_superfamily,
    molecular_family = Molecular_family
  )

### Long table of calls from every classifier version ---------------------------

bti <- read_tsv(file.path(input_dir, "bti_methylation_classifications.tsv"),
                show_col_types = FALSE) %>%
  rename(Kids_First_Biospecimen_ID = Bioassay_ID)

hist_methyl <- read_tsv(
  file.path(root_dir, "data", "histologies-base.tsv"),
  guess_max = 100000, show_col_types = FALSE
) %>%
  filter(experimental_strategy == "Methylation")

make_calls <- function(df, source, call, score, super_score = NA, class_score = NA) {
  tibble(
    Kids_First_Biospecimen_ID = df$Kids_First_Biospecimen_ID,
    source = source,
    call = df[[call]],
    score = df[[score]],
    super_score = if (is.na(super_score)) NA_real_ else df[[super_score]],
    class_score = if (is.na(class_score)) NA_real_ else df[[class_score]]
  )
}

all_calls <- bind_rows(
  make_calls(bti, "dkfz_v12.8", "dkfz_v12.8_subclass", "dkfz_v12.8_subclass_maxscore"),
  make_calls(bti, "nih_v3.1", "NIH_v3.1_subclass", "NIH_v3.1_subclass_meanscore",
             "NIH_v3.1_superfamily_meanscore", "NIH_v3.1_class_score"),
  make_calls(hist_methyl, "dkfz_v12", "dkfz_v12_methylation_subclass",
             "dkfz_v12_methylation_subclass_score"),
  make_calls(hist_methyl, "nih_v2", "NIH_v2_methylation_Class",
             "NIH_v2_methylation_Class_mean_score",
             "NIH_v2_methylation_Superfamily_mean_score",
             "NIH_v2_methylation_Class_mean_score")
) %>%
  filter(!is.na(call)) %>%
  mutate(
    key = norm_key(call),
    is_nih = grepl("^nih", source),
    # NIH calls only count as high confidence when the DKFZ call is not
    # (handled by precedence below), and both superfamily and class >= 0.9
    high_conf = if_else(
      is_nih,
      !is.na(super_score) & !is.na(class_score) &
        super_score >= nih_cutoff & class_score >= nih_cutoff,
      !is.na(score) & score >= dkfz_cutoff
    ),
    precedence = match(source, c("dkfz_v12.8", "nih_v3.1", "dkfz_v12", "nih_v2"))
  ) %>%
  left_join(key_table, by = "key") %>%
  # NIH control-tissue codes (CONTR_*, CTRL_*) are never subtyped, regardless of
  # what a crosswalk key resolves to (e.g., "CONTR HEMI" is listed as the NIH
  # v3 code for the GG row, which would otherwise call controls GNG, MAPK)
  mutate(is_control_code = grepl("^(CONTR|CTRL|NORM) ", key)) %>%
  left_join(crosswalk_info, by = "crosswalk_row") %>%
  mutate(
    no_call_class = no_call_class | is_control_code,
    OPC_molecular_subtype = if_else(is_control_code, NA_character_,
                                    OPC_molecular_subtype),
    rare = if_else(is_control_code, FALSE, rare)
  )

### Gaps: classes seen in data that the crosswalk can't resolve ------------------

gaps <- all_calls %>%
  filter(high_conf) %>%
  mutate(gap_reason = case_when(
    no_call_class %in% TRUE ~ NA_character_,
    is.na(crosswalk_row) ~ "class not in crosswalk",
    is.na(OPC_molecular_subtype) ~ "crosswalk row has no OPC_molecular_subtype"
  )) %>%
  filter(!is.na(gap_reason)) %>%
  count(source, call, gap_reason, name = "n_high_conf_biospecimens") %>%
  arrange(gap_reason, desc(n_high_conf_biospecimens))

write_tsv(gaps, file.path(results_dir, "methylation_crosswalk_gaps.tsv"))

### Resolve one call per biospecimen -------------------------------------------

# only calls that resolve to a crosswalk row can set the subtype; high-confidence
# calls that don't resolve are reported in the gaps table above
resolved <- all_calls %>%
  filter(high_conf, !is.na(crosswalk_row) | no_call_class) %>%
  # all no-call classes count as the same call when looking for conflicts
  mutate(call_identity = if_else(no_call_class, "no_call",
                                 as.character(crosswalk_row))) %>%
  group_by(Kids_First_Biospecimen_ID) %>%
  arrange(precedence, .by_group = TRUE) %>%
  mutate(
    # disagreement: another high-confidence call maps to a different crosswalk row
    n_distinct_rows = n_distinct(call_identity),
    conflicting_sources = paste(
      source[call_identity != first(call_identity)],
      collapse = ";"
    )
  ) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    Kids_First_Biospecimen_ID,
    methyl_confidence = "high",
    methyl_source = source,
    methyl_class = call,
    methyl_score = coalesce(class_score, score),
    OPC_molecular_subtype,
    methyl_no_call = no_call_class,
    rare,
    molecular_superfamily,
    molecular_family,
    methyl_conflict = n_distinct_rows > 1,
    methyl_conflicting_sources = if_else(n_distinct_rows > 1,
                                         conflicting_sources, NA_character_)
  )

# every methylation biospecimen gets a row, low / none confidence included
all_methyl_ids <- union(bti$Kids_First_Biospecimen_ID,
                        hist_methyl$Kids_First_Biospecimen_ID)

methyl_calls <- tibble(Kids_First_Biospecimen_ID = all_methyl_ids) %>%
  left_join(resolved, by = "Kids_First_Biospecimen_ID") %>%
  mutate(
    methyl_confidence = replace_na(methyl_confidence, "low_or_none"),
    methyl_conflict = replace_na(methyl_conflict, FALSE),
    methyl_no_call = replace_na(methyl_no_call, FALSE)
  ) %>%
  # high-confidence call exists in a class the crosswalk cannot resolve
  left_join(
    all_calls %>%
      filter(high_conf, is.na(crosswalk_row), !no_call_class) %>%
      distinct(Kids_First_Biospecimen_ID) %>%
      mutate(unmapped_high_conf_call = TRUE),
    by = "Kids_First_Biospecimen_ID"
  ) %>%
  mutate(unmapped_high_conf_call = replace_na(unmapped_high_conf_call, FALSE))

stopifnot(!anyDuplicated(methyl_calls$Kids_First_Biospecimen_ID))

write_tsv(methyl_calls, file.path(results_dir, "methylation_calls.tsv"))

# quick summary to the log
message("Biospecimens with methylation: ", nrow(methyl_calls))
message("High confidence: ", sum(methyl_calls$methyl_confidence == "high"))
message("High confidence, deliberate no-call (control / non-tumor): ",
        sum(methyl_calls$methyl_no_call))
message("High confidence, class has no OPC subtype yet: ",
        sum(methyl_calls$methyl_confidence == "high" &
              !methyl_calls$methyl_no_call &
              is.na(methyl_calls$OPC_molecular_subtype)))
message("Only unmappable high-confidence calls (see gaps table): ",
        sum(methyl_calls$methyl_confidence == "low_or_none" &
              methyl_calls$unmapped_high_conf_call))
message("Conflicting high-confidence calls across versions: ",
        sum(methyl_calls$methyl_conflict))
