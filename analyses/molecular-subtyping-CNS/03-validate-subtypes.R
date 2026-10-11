################################################################################
# 03-validate-subtypes.R
#
# Validate methylation-derived subtypes (results/methylation_calls.tsv) against
#   1. genetic lesions: input/subtype_lesion_rules.tsv, using SNVs (parquet from
#      02-prepare-lesion-tables.R) and fusions
#   2. historical pathology diagnosis: input/pathology_dx_groups.tsv and
#      input/subtype_family_dx_compatibility.tsv
#
# Methylation is never overridden. Lesions can only raise confidence, and a
# diagnosis that does not fit the subtype is flagged as a probable
# reclassification for review, since many historical diagnoses were wrong.
#
# Tumors with no high-confidence methylation call and a rule marked
# use_as_fallback are subtyped from lesions alone (lesion-only evidence).
#
# Usage: Rscript --vanilla 03-validate-subtypes.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS")
input_dir <- file.path(module_dir, "input")
results_dir <- file.path(module_dir, "results")
scratch_dir <- file.path(root_dir, "scratch", "molecular-subtyping-CNS")

### Inputs ---------------------------------------------------------------------

calls <- read_tsv(file.path(results_dir, "methylation_calls.tsv"),
                  show_col_types = FALSE)
rules <- read_tsv(file.path(input_dir, "subtype_lesion_rules.tsv"),
                  show_col_types = FALSE) %>%
  mutate(rule_order = row_number(),
         use_as_fallback = as.logical(use_as_fallback))
dx_groups <- read_tsv(file.path(input_dir, "pathology_dx_groups.tsv"),
                      show_col_types = FALSE)
family_compat <- read_tsv(file.path(input_dir, "subtype_family_dx_compatibility.tsv"),
                          show_col_types = FALSE) %>%
  arrange(order)

hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID, sample_id,
         match_id, cohort, sample_type, experimental_strategy, RNA_library,
         pathology_diagnosis)

snv_file <- file.path(scratch_dir, "snv_lesion_genes.parquet")
if (!file.exists(snv_file)) {
  stop("Run 02-prepare-lesion-tables.R first: ", snv_file, " not found")
}
con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
snv <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", snv_file)) %>%
  as_tibble() %>%
  # keep protein-altering variants only; hotspot rows without a protein change
  # can't match a pattern anyway
  filter(!is.na(HGVSp_Short)) %>%
  distinct()

fusion <- read_tsv(file.path(root_dir, "data", "fusion-annoFuse.tsv.gz"),
                   show_col_types = FALSE) %>%
  distinct(Kids_First_Biospecimen_ID = Sample, FusionName)

### Biospecimen groups by match_id ---------------------------------------------
# methylation is keyed on its own biospecimen; SNVs and fusions are keyed on the
# DNA and RNA biospecimens of the same tumor sample, linked through match_id

dna_strategies <- c("WGS", "WXS", "Targeted Sequencing")

dna_ids <- hist %>%
  filter(experimental_strategy %in% dna_strategies, !is.na(match_id)) %>%
  distinct(match_id, dna_id = Kids_First_Biospecimen_ID)
rna_ids <- hist %>%
  filter(experimental_strategy == "RNA-Seq", !is.na(match_id)) %>%
  distinct(match_id, rna_id = Kids_First_Biospecimen_ID)

### Lesion hits per rule -------------------------------------------------------

snv_rules <- rules %>% filter(lesion_type == "snv")
fusion_rules <- rules %>% filter(lesion_type == "fusion")
arm_rules <- rules %>% filter(lesion_type == "arm")

# GISTIC arm values (WGS only) for 1p/19q codeletion
arm_values <- read_tsv(file.path(scratch_dir, "arm_1p19q.tsv"),
                       show_col_types = FALSE)

snv_hits <- snv_rules %>%
  mutate(gene_regex = paste0("^(", gene, ")$")) %>%
  rowwise() %>%
  mutate(hits = list(
    snv %>%
      filter(str_detect(Hugo_Symbol, gene_regex),
             str_detect(HGVSp_Short, pattern)) %>%
      transmute(dna_id = Kids_First_Biospecimen_ID,
                evidence = paste0(Hugo_Symbol, " ", HGVSp_Short))
  )) %>%
  ungroup() %>%
  select(rule_id, OPC_molecular_subtype, hits) %>%
  unnest(hits) %>%
  inner_join(dna_ids, by = "dna_id", relationship = "many-to-many") %>%
  distinct(match_id, rule_id, OPC_molecular_subtype, evidence)

fusion_hits <- fusion_rules %>%
  rowwise() %>%
  mutate(hits = list(
    fusion %>%
      filter(str_detect(FusionName, pattern)) %>%
      transmute(rna_id = Kids_First_Biospecimen_ID,
                evidence = FusionName)
  )) %>%
  ungroup() %>%
  select(rule_id, OPC_molecular_subtype, hits) %>%
  unnest(hits) %>%
  inner_join(rna_ids, by = "rna_id", relationship = "many-to-many") %>%
  distinct(match_id, rule_id, OPC_molecular_subtype, evidence)

# an arm rule is met when every arm in `gene` (";" separated) is at or below the
# cutoff in `pattern`
arm_hits <- arm_rules %>%
  rowwise() %>%
  mutate(hits = list({
    arms <- str_split(gene, ";")[[1]]
    cutoff <- as.numeric(pattern)
    arm_values %>%
      filter(if_all(all_of(arms), ~ !is.na(.x) & .x <= cutoff)) %>%
      transmute(dna_id = Kids_First_Biospecimen_ID,
                evidence = paste0(paste(arms, collapse = "+"), " codeleted"))
  })) %>%
  ungroup() %>%
  select(rule_id, OPC_molecular_subtype, hits) %>%
  unnest(hits) %>%
  inner_join(dna_ids, by = "dna_id", relationship = "many-to-many") %>%
  distinct(match_id, rule_id, OPC_molecular_subtype, evidence)

lesion_hits <- bind_rows(snv_hits, fusion_hits, arm_hits)

# TP53 annotation (not part of the subtype), from tp53_nf1_score. Besides the
# module's tp53_altered status and classifier score, report the number of hits
# the way that module counts them: each SNV/indel is a hit, and a CNV loss or
# SV is one more hit (CNV and SV are treated as the same event). Two or more
# hits suggest both alleles are affected. One row per match_id, since DNA and RNA
# biospecimens of the same tumor each carry a row.
tp53 <- read_tsv(
  file.path(root_dir, "analyses", "tp53_nf1_score", "results",
            "tp53_altered_status.tsv"),
  col_types = cols(.default = col_character()), na = c("", "NA")
) %>%
  filter(!is.na(match_id)) %>%
  mutate(across(c(SNV_indel_counts, CNV_loss_counts, SV_counts, hotspot,
                  tp53_score), as.numeric),
         li_fraumeni = str_detect(replace_na(cancer_predispositions, ""),
                                  "Li-Fraumeni syndrome")) %>%
  group_by(match_id) %>%
  summarise(
    tp53_altered = case_when(
      any(tp53_altered == "loss") ~ "loss",
      any(tp53_altered == "activated") ~ "activated",
      any(tp53_altered == "Other") ~ "Other",
      TRUE ~ NA_character_
    ),
    tp53_classifier_score = suppressWarnings(max(tp53_score, na.rm = TRUE)),
    snv_hits = suppressWarnings(max(SNV_indel_counts, na.rm = TRUE)),
    cnv_sv_hit = suppressWarnings(max(CNV_loss_counts, SV_counts,
                                      na.rm = TRUE)) >= 1,
    tp53_hotspot = suppressWarnings(max(hotspot, na.rm = TRUE)) == 1,
    tp53_li_fraumeni = any(li_fraumeni),
    tp53_snv_evidence = paste(na.omit(unique(HGVSp_Short[HGVSp_Short != "NA"])),
                              collapse = "; "),
    tp53_cnv_evidence = paste(na.omit(unique(CNV_loss_evidence)),
                              collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    across(c(tp53_classifier_score, snv_hits), ~ if_else(is.infinite(.x), NA_real_, .x)),
    snv_hits = replace_na(snv_hits, 0),
    tp53_hotspot = replace_na(tp53_hotspot, FALSE),
    cnv_sv_hit = replace_na(cnv_sv_hit, FALSE),
    tp53_n_hits = snv_hits + as.integer(cnv_sv_hit),
    tp53_two_hits = tp53_n_hits >= 2
  ) %>%
  select(match_id, tp53_altered, tp53_classifier_score, tp53_n_hits,
         tp53_two_hits, tp53_hotspot, tp53_li_fraumeni, tp53_snv_evidence,
         tp53_cnv_evidence)

### Diagnosis helpers ----------------------------------------------------------

# subtype -> family and compatible dx groups (first match wins, top to bottom)
subtype_family <- function(subtype) {
  map_chr(subtype, function(s) {
    hit <- family_compat %>% filter(str_detect(s, subtype_regex)) %>% slice(1)
    if (nrow(hit) == 0) NA_character_ else hit$subtype_family
  })
}
family_dx <- function(family) {
  map_chr(family, function(f) {
    hit <- family_compat %>% filter(subtype_family == f) %>% slice(1)
    if (nrow(hit) == 0) NA_character_ else hit$compatible_dx_groups
  })
}

# a pathology_diagnosis can hold several values separated by ";"
dx_to_groups <- function(dx) {
  map_chr(dx, function(d) {
    if (is.na(d)) return(NA_character_)
    parts <- str_split(d, ";")[[1]] %>% str_trim()
    grp <- dx_groups$dx_group[match(parts, dx_groups$pathology_diagnosis)]
    paste(sort(unique(replace_na(grp, "UNMAPPED"))), collapse = ";")
  })
}

dx_concordance <- function(subtype, dx_group) {
  fam <- subtype_family(subtype)
  compat <- family_dx(fam)
  pmap_chr(list(fam, compat, dx_group), function(f, c, g) {
    if (is.na(g) || g %in% c("UNINFORMATIVE", "UNMAPPED")) return("no_informative_dx")
    if (is.na(f)) return("subtype_family_unmapped")
    allowed <- str_split(c, ";")[[1]]
    got <- str_split(g, ";")[[1]] %>% setdiff(c("UNINFORMATIVE", "UNMAPPED"))
    if (length(got) == 0) return("no_informative_dx")
    primary <- allowed[1]
    if (primary %in% got) "concordant"
    else if (any(got %in% allowed)) "compatible"
    else "discordant"
  })
}

### Methylation-based subtypes: validate ---------------------------------------

methyl <- calls %>%
  filter(methyl_confidence == "high", !methyl_control,
         !is.na(OPC_molecular_subtype)) %>%
  left_join(hist %>% select(Kids_First_Biospecimen_ID,
                            Kids_First_Participant_ID, sample_id, match_id,
                            pathology_diagnosis),
            by = "Kids_First_Biospecimen_ID")

# lesions confirming this subtype (rules for the same subtype only)
confirming <- lesion_hits %>%
  group_by(match_id, OPC_molecular_subtype) %>%
  summarise(lesion_evidence = paste(sort(unique(evidence)), collapse = "; "),
            .groups = "drop")

# does this tumor have the assay a rule needs?
has_dna <- unique(dna_ids$match_id)
has_rna <- unique(rna_ids$match_id)
rule_assays <- rules %>%
  group_by(OPC_molecular_subtype) %>%
  summarise(needs_dna = any(lesion_type %in% c("snv", "arm")),
            needs_rna = any(lesion_type == "fusion"),
            has_rule = TRUE, .groups = "drop")

validated <- methyl %>%
  left_join(rule_assays, by = "OPC_molecular_subtype") %>%
  left_join(confirming, by = c("match_id", "OPC_molecular_subtype")) %>%
  mutate(
    assay_available = (needs_dna %in% TRUE & match_id %in% has_dna) |
      (needs_rna %in% TRUE & match_id %in% has_rna),
    lesion_status = case_when(
      !(has_rule %in% TRUE) ~ "no_rule",
      !is.na(lesion_evidence) ~ "confirmed",
      assay_available ~ "not_detected",
      TRUE ~ "no_assay"
    ),
    dx_group = dx_to_groups(pathology_diagnosis),
    dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group),
    evidence_level = case_when(
      lesion_status == "confirmed" ~ "methylation + lesion",
      TRUE ~ "methylation only"
    ),
    needs_review = dx_concordance == "discordant" | methyl_conflict
  ) %>%
  left_join(tp53, by = "match_id") %>%
  select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID, sample_id,
         match_id, OPC_molecular_subtype, methyl_class, methyl_source,
         methyl_score, methyl_conflict, methyl_conflicting_sources,
         lesion_status, lesion_evidence, pathology_diagnosis, dx_group,
         dx_concordance, tp53_altered, tp53_classifier_score, tp53_n_hits,
         tp53_two_hits, tp53_hotspot, tp53_li_fraumeni, tp53_snv_evidence,
         tp53_cnv_evidence, evidence_level, needs_review)

write_tsv(validated, file.path(results_dir, "methylation_subtype_validation.tsv"))
validated %>%
  filter(needs_review) %>%
  arrange(dx_concordance, OPC_molecular_subtype) %>%
  write_tsv(file.path(results_dir, "subtypes_for_review.tsv"))

### Lesion-only fallback for tumors without a high-confidence methylation call -

methyl_match_ids <- validated$match_id %>% unique()
control_match_ids <- calls %>%
  filter(methyl_confidence == "high", methyl_control) %>%
  inner_join(hist %>% select(Kids_First_Biospecimen_ID, match_id),
             by = "Kids_First_Biospecimen_ID") %>%
  pull(match_id) %>% unique()

fallback_rule_ids <- rules %>% filter(use_as_fallback) %>% pull(rule_id)

# fallback only applies to CNS tumors (PBTA cohort); e.g., IDH1 R132 in TARGET
# leukemia is not an HGG, IDH call. Tumors whose methylation is a control /
# low-purity tissue are excluded because their lesions are not reliable.
pbta_match_ids <- hist %>% filter(cohort == "PBTA") %>% pull(match_id) %>% unique()

fallback <- lesion_hits %>%
  filter(rule_id %in% fallback_rule_ids,
         match_id %in% pbta_match_ids,
         !match_id %in% methyl_match_ids,
         !match_id %in% control_match_ids) %>%
  left_join(rules %>% select(rule_id, rule_order), by = "rule_id") %>%
  group_by(match_id) %>%
  arrange(rule_order, .by_group = TRUE) %>%
  summarise(OPC_molecular_subtype = first(OPC_molecular_subtype),
            lesion_evidence = paste(unique(evidence), collapse = "; "),
            n_distinct_subtypes = n_distinct(OPC_molecular_subtype),
            .groups = "drop") %>%
  left_join(hist %>%
              filter(sample_type == "Tumor") %>%
              distinct(match_id, Kids_First_Participant_ID, sample_id,
                       pathology_diagnosis),
            by = "match_id") %>%
  mutate(dx_group = dx_to_groups(pathology_diagnosis))

# IDH-mutant gliomas: 1p/19q codeletion makes an IDH-mutant oligodendroglioma.
# Without codeletion evidence the grade comes from the diagnosis group, and an
# oligodendroglioma diagnosis with no codeletion evidence is held unclassified.
codeleted_match_ids <- arm_hits %>% filter(rule_id == "OLIGO_1P19Q") %>%
  pull(match_id) %>% unique()
# tumors with arm-level data at all (WGS GISTIC), so "not codeleted" can be told
# apart from "not assessable"
arm_assessed_match_ids <- dna_ids %>%
  filter(dna_id %in% arm_values$Kids_First_Biospecimen_ID) %>%
  pull(match_id) %>% unique()

fallback <- fallback %>%
  mutate(
    is_idh = OPC_molecular_subtype == "HGG, IDH",
    codeleted = match_id %in% codeleted_match_ids,
    oligo_dx = str_detect(replace_na(pathology_diagnosis, ""),
                          regex("oligodendroglioma", ignore_case = TRUE)),
    hold_reason = case_when(
      is_idh & !codeleted & oligo_dx & match_id %in% arm_assessed_match_ids ~
        "oligodendroglioma diagnosis, 1p/19q not codeleted",
      is_idh & !codeleted & oligo_dx ~
        "oligodendroglioma diagnosis, 1p/19q not assessable (no WGS arm data)",
      TRUE ~ NA_character_
    ),
    OPC_molecular_subtype = case_when(
      is_idh & codeleted ~ "Oligo, IDH",
      is_idh & !codeleted & !oligo_dx & str_detect(dx_group, "^LGG$") ~ "LGG, IDH",
      TRUE ~ OPC_molecular_subtype
    ),
    lesion_evidence = if_else(is_idh & codeleted,
                              paste0(lesion_evidence, "; 1p+19q codeleted"),
                              lesion_evidence)
  ) %>%
  select(-is_idh, -oligo_dx) %>%
  left_join(tp53, by = "match_id") %>%
  mutate(
    dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group),
    evidence_level = "lesion only",
    # a lesion-only call that the diagnosis contradicts, or a held call, is not
    # assigned
    assign = dx_concordance != "discordant" & is.na(hold_reason)
  )

write_tsv(fallback, file.path(results_dir, "lesion_fallback_subtypes.tsv"))

# pathology diagnoses the dx group table does not recognize
unmapped_fb_dx <- fallback %>%
  filter(str_detect(replace_na(dx_group, ""), "UNMAPPED")) %>%
  distinct(pathology_diagnosis)

### Summary --------------------------------------------------------------------

message("\n== Methylation-based subtypes: ", nrow(validated))
print(validated %>% count(lesion_status))
print(validated %>% count(dx_concordance))
message("Flagged for review: ", sum(validated$needs_review))

message("\n== Lesion-only fallback: ", nrow(fallback), " tumors without methylation")
print(fallback %>% count(OPC_molecular_subtype, assign, hold_reason))

# subtype families the diagnosis table does not recognize
unmapped_family <- validated %>%
  filter(dx_concordance == "subtype_family_unmapped") %>%
  distinct(OPC_molecular_subtype)
if (nrow(unmapped_family) > 0) {
  message("\nSubtypes with no family in subtype_family_dx_compatibility.tsv:")
  print(unmapped_family$OPC_molecular_subtype)
}
unmapped_dx <- validated %>%
  filter(str_detect(replace_na(dx_group, ""), "UNMAPPED")) %>%
  distinct(pathology_diagnosis) %>%
  bind_rows(unmapped_fb_dx) %>%
  distinct()
if (nrow(unmapped_dx) > 0) {
  message("\nPathology diagnoses not in pathology_dx_groups.tsv:")
  print(unmapped_dx$pathology_diagnosis)
}
