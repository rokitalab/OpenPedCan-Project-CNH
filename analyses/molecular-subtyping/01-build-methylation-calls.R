################################################################################
# 01-build-methylation-calls.R
#
# Build one resolved, per-biospecimen methylation call table that downstream
# molecular subtyping modules consume instead of re-implementing thresholds and
# classifier-version aliasing.
#
# Calls come from input/bti_methylation_classifications.tsv: DKFZ v12.8 and
# NIH v3.1. Earlier classifier versions (DKFZ v11 / v12, NIH v2) are not used.
# The serum-free derived cell lines have no array in the release; their classes are
# those of the CHOP DGD classifier (input/serum_free_cell_line_methylation.tsv, one
# row per sample), attached to the serum-free cell lines of the sample.
#
# Confidence rule (same as molecular-subtyping-Rare-CNS):
#   high  = DKFZ subclass score >= 0.8, or (DKFZ missing / < 0.8) NIH
#           superfamily and class scores >= 0.9
#   CHOP DGD (cell lines): class mean score >= 0.8
# Precedence among passing calls: DKFZ v12.8, then NIH v3.1, then CHOP DGD
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
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
input_dir <- file.path(module_dir, "input")
results_dir <- file.path(module_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

dkfz_cutoff <- 0.8
nih_cutoff <- 0.9
chop_dgd_cutoff <- 0.8

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

# Subtype names for classes that have no OPC_molecular_subtype in the crosswalk
# yet, kept as a reviewable table with the rationale. Existing crosswalk values
# are never overwritten.
additions <- read_tsv(file.path(input_dir, "opc_subtype_additions.tsv"),
                      show_col_types = FALSE, na = character()) %>%
  select(Abbreviation_internal_NIHv2, subtype_addition = OPC_molecular_subtype,
         non_cns)
stopifnot(!anyDuplicated(additions$Abbreviation_internal_NIHv2),
          all(additions$Abbreviation_internal_NIHv2 %in%
                crosswalk$Abbreviation_internal_NIHv2))

crosswalk <- crosswalk %>%
  left_join(additions, by = "Abbreviation_internal_NIHv2") %>%
  mutate(
    subtype_addition = if_else(!is.na(OPC_molecular_subtype), NA_character_,
                               subtype_addition),
    OPC_molecular_subtype = coalesce(OPC_molecular_subtype, subtype_addition),
    # classes that are not primary CNS tumors (kept so every classifier call
    # resolves, but flagged for downstream exclusion or review)
    non_cns = replace_na(as.logical(non_cns), FALSE)
  )

key_cols <- c("dkfz_Abbreviation", "dkfz_Abbreviation_v12.8",
              "Abbreviation_internal_NIHv2", "Abbreviation_internal_NIHv3",
              "OPC_molecular_subtype_methyl")

key_table <- crosswalk %>%
  select(crosswalk_row, all_of(key_cols)) %>%
  pivot_longer(-crosswalk_row, values_to = "raw_key") %>%
  mutate(key = norm_key(raw_key)) %>%
  filter(!is.na(key)) %>%
  distinct(key, crosswalk_row)

# NIH v3.1 subclass names that the crosswalk does not list, mapped to the
# crosswalk row of the DKFZ v12.8 subclass they match one to one in the arrays
# where both classifiers are high confidence (input/nih_v3.1_class_aliases.tsv)
nih_aliases <- read_tsv(file.path(input_dir, "nih_v3.1_class_aliases.tsv"),
                        show_col_types = FALSE) %>%
  mutate(alias_key = norm_key(`nih_v3.1_subclass`),
         dkfz_key = norm_key(`dkfz_v12.8_subclass`))
alias_keys <- nih_aliases %>%
  left_join(key_table %>% rename(dkfz_key = key), by = "dkfz_key",
            relationship = "many-to-many")
stopifnot(!anyNA(alias_keys$crosswalk_row))
key_table <- bind_rows(
  key_table,
  alias_keys %>% transmute(key = alias_key, crosswalk_row)
) %>%
  distinct(key, crosswalk_row)

# classes seen in the classifier output that the crosswalk does not list yet get
# a synthetic crosswalk row (input/extra_methylation_classes.tsv)
extra_classes <- read_tsv(file.path(input_dir, "extra_methylation_classes.tsv"),
                          show_col_types = FALSE) %>%
  mutate(crosswalk_row = 10000L + row_number())

key_table <- bind_rows(
  key_table,
  extra_classes %>% transmute(key = norm_key(class_key), crosswalk_row)
)

# a key must resolve to exactly one crosswalk row, otherwise calls are ambiguous
ambiguous_keys <- key_table %>%
  count(key) %>%
  filter(n > 1)
if (nrow(ambiguous_keys) > 0) {
  warning("Ambiguous crosswalk keys dropped: ",
          paste(ambiguous_keys$key, collapse = ", "))
  key_table <- key_table %>% filter(!key %in% ambiguous_keys$key)
}

# the crosswalk marks control classes (control tissues, inflammatory
# microenvironment; low tumor purity) with "Control" in OPC_molecular_subtype;
# a blank cell instead means the class has not been assigned a subtype yet
crosswalk_info <- bind_rows(
  crosswalk %>% select(crosswalk_row, OPC_molecular_subtype, Rare,
                       Molecular_superfamily, Molecular_family, non_cns),
  extra_classes %>% transmute(crosswalk_row, OPC_molecular_subtype,
                              Rare = NA_character_, Molecular_superfamily = NA_character_,
                              Molecular_family = NA_character_, non_cns = FALSE)
) %>%
  transmute(
    crosswalk_row,
    no_call_class = OPC_molecular_subtype %in% c("Control", "NA"),
    OPC_molecular_subtype = if_else(no_call_class, NA_character_,
                                    OPC_molecular_subtype),
    rare = Rare %in% "Yes",
    non_cns,
    molecular_superfamily = Molecular_superfamily,
    molecular_family = Molecular_family
  )

### Long table of calls from every classifier version ---------------------------

bti <- read_tsv(file.path(input_dir, "bti_methylation_classifications.tsv"),
                show_col_types = FALSE) %>%
  rename(Kids_First_Biospecimen_ID = Bioassay_ID)

# every methylation biospecimen in the histologies (ids only: the calls are not
# read from there), so arrays without a v12.8 / v3.1 classification are listed
# as low_or_none
hist_methyl <- read_tsv(
  file.path(root_dir, "data", "histologies-base.tsv"),
  col_select = c(Kids_First_Biospecimen_ID, experimental_strategy),
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

### CHOP DGD classes of the serum-free cell lines ------------------------------

# The file has one row per sample (no passage), so the class is attached to every
# serum-free cell line of the sample (sample_id and a serum-free composition). The
# lines have no array, so the call is set on one biospecimen of each line (the WGS /
# WXS one when there is one): downstream steps read a call by its biospecimen and
# take the line's own match_id, DNA and RNA from the histologies.
chop_dgd <- read_tsv(file.path(input_dir, "serum_free_cell_line_methylation.tsv"),
                     show_col_types = FALSE, col_types = cols(.default = "c")) %>%
  mutate(class_code = sub("^Methylation_class_", "", class),
         across(c(class_mean_score, class_consistency_score), as.numeric))
stopifnot(!anyDuplicated(chop_dgd$sample_id))
# CHOP DGD codes that the crosswalk spells differently (input/chop_dgd_class_aliases.tsv)
chop_dgd_aliases <- read_tsv(file.path(input_dir, "chop_dgd_class_aliases.tsv"),
                             show_col_types = FALSE)
chop_dgd <- chop_dgd %>%
  mutate(call_class = coalesce(
    chop_dgd_aliases$crosswalk_class[match(class_code, chop_dgd_aliases$chop_dgd_class)],
    class_code))

serum_free_lines <- read_tsv(
  file.path(root_dir, "data", "histologies-base.tsv"),
  col_select = c(Kids_First_Biospecimen_ID, sample_id, match_id, composition,
                 cell_line_composition, cell_line_passage, experimental_strategy),
  guess_max = 100000, show_col_types = FALSE) %>%
  filter(composition == "Derived Cell Line",
         str_starts(replace_na(cell_line_composition, ""), "Serum-free"),
         sample_id %in% chop_dgd$sample_id) %>%
  mutate(assay_order = match(experimental_strategy, c("WGS", "WXS"), nomatch = 3L)) %>%
  arrange(match_id, assay_order, Kids_First_Biospecimen_ID) %>%
  group_by(match_id) %>%
  slice(1) %>%
  ungroup()

chop_dgd_lines <- chop_dgd %>%
  inner_join(serum_free_lines %>%
               select(sample_id, match_id, Kids_First_Biospecimen_ID,
                      cell_line_composition, cell_line_passage),
             by = "sample_id")
stopifnot(!anyDuplicated(chop_dgd_lines$Kids_First_Biospecimen_ID))

all_calls <- bind_rows(
  make_calls(bti, "dkfz_v12.8", "dkfz_v12.8_subclass", "dkfz_v12.8_subclass_maxscore"),
  make_calls(bti, "nih_v3.1", "NIH_v3.1_subclass", "NIH_v3.1_subclass_meanscore",
             "NIH_v3.1_superfamily_meanscore", "NIH_v3.1_class_score"),
  make_calls(chop_dgd_lines %>% mutate(call = call_class, score = class_mean_score),
             "chop_dgd", "call", "score")
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
      if_else(source == "chop_dgd", !is.na(score) & score >= chop_dgd_cutoff,
              !is.na(score) & score >= dkfz_cutoff)
    ),
    precedence = match(source, c("dkfz_v12.8", "nih_v3.1", "chop_dgd"))
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
  # all control classes count as the same call when looking for conflicts
  mutate(call_identity = if_else(no_call_class, "control",
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
    methyl_control = no_call_class,
    methyl_non_cns = non_cns,
    rare,
    molecular_superfamily,
    molecular_family,
    methyl_conflict = n_distinct_rows > 1,
    methyl_conflicting_sources = if_else(n_distinct_rows > 1,
                                         conflicting_sources, NA_character_)
  )

# every methylation biospecimen gets a row, low / none confidence included
all_methyl_ids <- union(union(bti$Kids_First_Biospecimen_ID,
                              hist_methyl$Kids_First_Biospecimen_ID),
                        chop_dgd_lines$Kids_First_Biospecimen_ID)

methyl_calls <- tibble(Kids_First_Biospecimen_ID = all_methyl_ids) %>%
  left_join(resolved, by = "Kids_First_Biospecimen_ID") %>%
  mutate(
    methyl_confidence = replace_na(methyl_confidence, "low_or_none"),
    methyl_conflict = replace_na(methyl_conflict, FALSE),
    methyl_control = replace_na(methyl_control, FALSE),
    methyl_non_cns = replace_na(methyl_non_cns, FALSE)
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

# the CHOP DGD classes and where they went: one row per sample and serum-free cell
# line; a sample with no serum-free line in the histologies has no match_id
chop_dgd_matching <- chop_dgd %>%
  left_join(chop_dgd_lines %>%
              select(sample_id, match_id, Kids_First_Biospecimen_ID,
                     cell_line_composition, cell_line_passage),
            by = "sample_id") %>%
  left_join(methyl_calls %>%
              select(Kids_First_Biospecimen_ID, methyl_confidence, methyl_class,
                     OPC_molecular_subtype, methyl_control),
            by = "Kids_First_Biospecimen_ID") %>%
  transmute(sample_id, match_id, Kids_First_Biospecimen_ID, cell_line_composition,
            cell_line_passage, chop_dgd_class = class_code,
            class_mean_score, class_consistency_score,
            result = case_when(
              is.na(match_id) ~ "no serum-free cell line in the histologies",
              methyl_confidence != "high" ~ "class mean score below 0.8",
              methyl_control ~ "control class",
              is.na(OPC_molecular_subtype) ~ "class not in the crosswalk or without subtype",
              TRUE ~ "subtype from methylation"),
            methyl_subtype = OPC_molecular_subtype)
write_tsv(chop_dgd_matching, file.path(results_dir, "chop_dgd_cell_line_calls.tsv"))

# quick summary to the log
message("Biospecimens with methylation: ", nrow(methyl_calls))
message("High confidence: ", sum(methyl_calls$methyl_confidence == "high"))
message("High confidence, control (low tumor purity): ",
        sum(methyl_calls$methyl_control))
message("High confidence, class has no OPC subtype yet: ",
        sum(methyl_calls$methyl_confidence == "high" &
              !methyl_calls$methyl_control &
              is.na(methyl_calls$OPC_molecular_subtype)))
message("Only unmappable high-confidence calls (see gaps table): ",
        sum(methyl_calls$methyl_confidence == "low_or_none" &
              methyl_calls$unmapped_high_conf_call))
message("Conflicting high-confidence calls across versions: ",
        sum(methyl_calls$methyl_conflict))
