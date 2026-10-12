################################################################################
# 03-validate-subtypes.R
#
# Validate methylation-derived subtypes (results/methylation_calls.tsv) against
#   1. genetic lesions (input/subtype_lesion_rules.tsv)
#   2. the historical pathology diagnosis (input/pathology_dx_groups.tsv and
#      input/subtype_family_dx_compatibility.tsv)
#
# Methylation is never overridden: lesions raise confidence, and a diagnosis
# that does not fit the subtype is flagged as a probable reclassification.
# PBTA tumors with no high-confidence methylation call are subtyped from
# lesions alone for rules marked use_as_fallback.
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

source(file.path(module_dir, "utils", "lesion-hits.R"))
source(file.path(module_dir, "utils", "dx-concordance.R"))
source(file.path(module_dir, "utils", "tp53-annotation.R"))
source(file.path(module_dir, "utils", "gene-annotation.R"))

### Inputs ---------------------------------------------------------------------

calls <- read_tsv(file.path(results_dir, "methylation_calls.tsv"),
                  show_col_types = FALSE)
rules <- read_tsv(file.path(input_dir, "subtype_lesion_rules.tsv"),
                  show_col_types = FALSE) %>%
  mutate(rule_order = row_number(),
         use_as_fallback = as.logical(use_as_fallback),
         fallback_requires_dx = as.logical(fallback_requires_dx),
         definitive = as.logical(definitive),
         lesion_overrides_methylation = as.logical(lesion_overrides_methylation))
dx_groups <- read_tsv(file.path(input_dir, "pathology_dx_groups.tsv"),
                      show_col_types = FALSE)
family_compat <- read_tsv(
  file.path(input_dir, "subtype_family_dx_compatibility.tsv"),
  show_col_types = FALSE
) %>%
  arrange(order)

hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID, sample_id,
         match_id, cohort, sample_type, experimental_strategy,
         pathology_diagnosis, pathology_free_text_diagnosis)

snv_file <- file.path(scratch_dir, "snv_lesion_genes.parquet")
if (!file.exists(snv_file)) {
  stop("Run 02-prepare-lesion-tables.R first: ", snv_file, " not found")
}
con <- dbConnect(duckdb::duckdb())
snv <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", snv_file)) %>%
  as_tibble() %>%
  filter(!is.na(HGVSp_Short)) %>%
  distinct()
dbDisconnect(con, shutdown = TRUE)

fusion <- read_tsv(file.path(root_dir, "data", "fusion-annoFuse.tsv.gz"),
                   show_col_types = FALSE) %>%
  distinct(Kids_First_Biospecimen_ID = Sample, FusionName)

arm_values <- read_tsv(file.path(scratch_dir, "arm_1p19q.tsv"),
                       show_col_types = FALSE)

con <- dbConnect(duckdb::duckdb())
cnv <- dbGetQuery(con, sprintf(
  "SELECT * FROM read_parquet('%s')",
  file.path(scratch_dir, "cnv_lesion_genes.parquet")
)) %>% as_tibble()
dbDisconnect(con, shutdown = TRUE)

tp53 <- read_tp53_annotation(root_dir)
annotation_genes <- read_tsv(file.path(input_dir, "annotation_genes.tsv"),
                             show_col_types = FALSE)
expression_annotations <- read_tsv(
  file.path(input_dir, "expression_annotations.tsv"), show_col_types = FALSE
)
expression_long <- read_tsv(
  file.path(scratch_dir, "expression_annotation_genes.tsv"),
  show_col_types = FALSE
)
cnv_segments <- read_tsv(
  file.path(root_dir, "data", "consensus_seg_with_status.tsv"),
  col_types = cols_only(Kids_First_Biospecimen_ID = "c", chrom = "c",
                        loc.start = "d", loc.end = "d", copy.num = "d",
                        tumor_ploidy = "d")
)

### Lesion hits ----------------------------------------------------------------
# methylation is keyed on its own biospecimen, while SNVs, arm values and
# fusions are keyed on the DNA / RNA biospecimens of the same tumor, so they
# are linked through match_id

assay_ids <- function(strategies) {
  hist %>%
    filter(experimental_strategy %in% strategies, !is.na(match_id)) %>%
    distinct(match_id, id = Kids_First_Biospecimen_ID)
}
dna_ids <- assay_ids(c("WGS", "WXS", "Targeted Sequencing"))
rna_ids <- assay_ids("RNA-Seq")

lesion_hits <- find_lesion_hits(
  rules, snv, fusion, arm_values, cnv,
  assay_ids = list(snv = dna_ids, arm = dna_ids, cnv = dna_ids, fusion = rna_ids)
)

# H3 K27M/I and G34/G35 R/V hits, for the h3_status annotation (H3 wildtype
# means a DNA assay with none of these)
h3_hits <- lesion_hits %>%
  filter(rule_id %in% c("DMG_H3K28", "DHG_H3G35", "DMG_EGFR_H3K28"))

### Methylation-based subtypes -------------------------------------------------

rule_assays <- rules %>%
  group_by(OPC_molecular_subtype) %>%
  summarise(needs_dna = any(lesion_type %in% c("snv", "arm", "cnv")),
            needs_rna = any(lesion_type == "fusion"),
            has_rule = TRUE, .groups = "drop")

confirming <- lesion_hits %>%
  group_by(match_id, OPC_molecular_subtype) %>%
  summarise(lesion_evidence = paste(sort(unique(evidence)), collapse = "; "),
            .groups = "drop")

validated <- calls %>%
  filter(methyl_confidence == "high", !methyl_control,
         !is.na(OPC_molecular_subtype)) %>%
  left_join(hist %>% select(Kids_First_Biospecimen_ID,
                            Kids_First_Participant_ID, sample_id, match_id,
                            pathology_diagnosis, pathology_free_text_diagnosis),
            by = "Kids_First_Biospecimen_ID") %>%
  left_join(rule_assays, by = "OPC_molecular_subtype") %>%
  left_join(confirming, by = c("match_id", "OPC_molecular_subtype")) %>%
  mutate(
    assay_available = (needs_dna %in% TRUE & match_id %in% dna_ids$match_id) |
      (needs_rna %in% TRUE & match_id %in% rna_ids$match_id),
    lesion_status = case_when(
      !(has_rule %in% TRUE) ~ "no_rule",
      !is.na(lesion_evidence) ~ "confirmed",
      assay_available ~ "not_detected",
      TRUE ~ "no_assay"
    ),
    dx_group = dx_to_groups(pathology_diagnosis, dx_groups),
    # PXA (issue 113) is split by BRAF V600E (the PXA rule's lesion), for DKFZ
    # (score >= 0.8) and NIH (>= 0.9) high-confidence calls alike. Without a DNA
    # assay the status is not assessable, so the tumor stays plain PXA. The
    # HGG / LGG prefix follows the pathology diagnosis, or the free-text
    # diagnosis when that is missing or uninformative; LGG is the default
    # (issue 113).
    dx_informative = !is.na(dx_group) & !dx_group %in% c("UNINFORMATIVE", "UNMAPPED"),
    pxa_prefix = case_when(
      str_detect(replace_na(dx_group, ""), "HGG") ~ "HGG",
      dx_informative ~ "LGG",
      free_text_grade_group(pathology_free_text_diagnosis) %in% "HGG" ~ "HGG",
      TRUE ~ "LGG"
    ),
    OPC_molecular_subtype = case_when(
      OPC_molecular_subtype == "PXA" & lesion_status == "confirmed" ~
        paste0(pxa_prefix, ", PXA BRAF V600E"),
      OPC_molecular_subtype == "PXA" & lesion_status == "not_detected" ~
        paste0(pxa_prefix, ", PXA BRAF V600E negative"),
      OPC_molecular_subtype == "PXA" ~ paste0(pxa_prefix, ", PXA"),
      TRUE ~ OPC_molecular_subtype
    ),
    # "-mutant" is only used when the H3 variant was found. A methylation class
    # without that evidence keeps the plain label ("DMG, H3 K28", "DHG, H3 G35"),
    # and the DMG_EGFR class is H3 K28-altered, EGFR (or H3 K28-mutant, EGFR
    # with the variant).
    OPC_molecular_subtype = case_when(
      OPC_molecular_subtype == "DMG, H3 K28" & lesion_status == "confirmed" ~
        "DMG, H3 K28-mutant",
      OPC_molecular_subtype == "DHG, H3 G35" & lesion_status == "confirmed" ~
        "DHG, H3 G35-mutant",
      OPC_molecular_subtype == "DMG, H3 K28-altered, EGFR" &
        lesion_status == "confirmed" ~ "DMG, H3 K28-mutant, EGFR",
      TRUE ~ OPC_molecular_subtype
    ),
    dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                    family_compat),
    evidence_level = if_else(lesion_status == "confirmed",
                             "methylation + lesion", "methylation only"),
    needs_review = dx_concordance == "discordant" | methyl_conflict
  ) %>%
  left_join(tp53, by = "match_id") %>%
  annotate_genes(annotation_genes, snv, dna_ids, cnv_segments) %>%
  annotate_expression(expression_annotations, expression_long, rna_ids) %>%
  add_h3_status(h3_hits, dna_ids) %>%
  relabel_egfr_only() %>%
  add_h3k27_mechanism(h3_mutated = .$lesion_status == "confirmed",
                      has_dna = .$match_id %in% dna_ids$match_id,
                      has_rna = .$match_id %in% rna_ids$match_id) %>%
  select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID, sample_id,
         match_id, OPC_molecular_subtype, methyl_class, methyl_source,
         methyl_score, methyl_conflict, methyl_conflicting_sources,
         lesion_status, lesion_evidence, pathology_diagnosis, dx_group,
         dx_concordance, starts_with("tp53_"),
         ends_with("_mutation_status"), ends_with("_variants"),
         ends_with("_amplification"), ends_with("_activating_alteration"),
         ends_with("_tpm"), ends_with("_z"), ends_with("_expression_status"),
         starts_with("h3_"), starts_with("h3k27_"), evidence_level, needs_review)

# an H3 K27M/I or G34/G35 R/V variant defines DMG/DHG, so a methylation subtype
# that is not DMG/DHG contradicts it and is flagged for review
# A definitive lesion (CTNNB1 exon 3 for MB WNT / CRANIO ADAM, BRAF V600E for
# CRANIO PAP) that belongs to a different subtype of the same family
# contradicts the methylation call. In craniopharyngioma the lesion decides:
# CTNNB1 exon 3 and BRAF V600E are mutually exclusive, and pathology confirmed
# adamantinomatous in 5 of 5 reviewed methylation-PAP tumors (see README). For
# other families the methylation call is kept and the tumor flagged, unless the
# call is also weak (NIH-only, or a DKFZ score within 0.05 of the 0.8 cutoff),
# in which case it is set aside and the tumor is subtyped by the lesion
# fallback.
borderline_dkfz <- 0.85

definitive_hits <- lesion_hits %>%
  inner_join(select(rules, rule_id, definitive, lesion_overrides_methylation),
             by = "rule_id") %>%
  filter(definitive) %>%
  group_by(match_id, OPC_molecular_subtype, lesion_overrides_methylation) %>%
  summarise(evidence = paste(sort(unique(evidence)), collapse = "; "),
            .groups = "drop") %>%
  transmute(match_id, other_subtype = OPC_molecular_subtype,
            other_family = str_extract(OPC_molecular_subtype, "^[^,]+"),
            other_evidence = evidence, overrides = lesion_overrides_methylation)

contradictions <- validated %>%
  mutate(family = str_extract(OPC_molecular_subtype, "^[^,]+")) %>%
  inner_join(definitive_hits, by = "match_id", relationship = "many-to-many") %>%
  filter(family == other_family, OPC_molecular_subtype != other_subtype) %>%
  group_by(match_id) %>%
  summarise(
    lesion_contradiction = paste0(other_evidence, " (", other_subtype, ")",
                                  collapse = "; "),
    override_subtype = first(other_subtype[overrides]),
    override_evidence = first(other_evidence[overrides]),
    .groups = "drop"
  )

validated <- validated %>%
  left_join(contradictions, by = "match_id") %>%
  mutate(
    h3_variant_conflict = match_id %in% h3_hits$match_id &
      !str_detect(OPC_molecular_subtype, "^(DMG|DHG)"),
    weak_methylation = str_detect(methyl_source, "^nih") |
      methyl_score < borderline_dkfz,
    # the lesion decides: relabel, keep the methylation subtype for the record
    methylation_overridden = !is.na(override_subtype),
    methylation_subtype = OPC_molecular_subtype,
    OPC_molecular_subtype = coalesce(override_subtype, OPC_molecular_subtype),
    lesion_status = if_else(methylation_overridden, "confirmed", lesion_status),
    lesion_evidence = if_else(methylation_overridden, override_evidence,
                              lesion_evidence),
    evidence_level = if_else(methylation_overridden,
                             "lesion overrides methylation", evidence_level),
    # resolved by the rule, so no longer a contradiction to review
    lesion_contradiction = if_else(methylation_overridden, NA_character_,
                                   lesion_contradiction),
    methylation_downgraded = !is.na(lesion_contradiction) & weak_methylation,
    needs_review = if_else(
      methylation_overridden,
      dx_concordance == "discordant" | methyl_conflict | h3_variant_conflict,
      needs_review | h3_variant_conflict | !is.na(lesion_contradiction)
    )
  )

validated %>%
  filter(methylation_downgraded) %>%
  select(Kids_First_Biospecimen_ID, sample_id, match_id, methyl_source,
         methyl_class, methyl_score, methylation_subtype = OPC_molecular_subtype,
         lesion_contradiction, pathology_diagnosis) %>%
  write_tsv(file.path(results_dir, "methylation_downgraded.tsv"))

validated %>%
  filter(methylation_overridden) %>%
  select(Kids_First_Biospecimen_ID, sample_id, match_id, methyl_source,
         methyl_class, methyl_score, methylation_subtype, lesion_subtype = OPC_molecular_subtype,
         lesion_evidence, pathology_diagnosis) %>%
  write_tsv(file.path(results_dir, "methylation_overridden_by_lesion.tsv"))

validated <- validated %>%
  filter(!methylation_downgraded) %>%
  select(-weak_methylation, -methylation_downgraded, -override_subtype,
         -override_evidence)

write_tsv(validated, file.path(results_dir, "methylation_subtype_validation.tsv"))
validated %>%
  filter(needs_review) %>%
  arrange(dx_concordance, OPC_molecular_subtype) %>%
  write_tsv(file.path(results_dir, "subtypes_for_review.tsv"))

### Lesion-only fallback -------------------------------------------------------
# PBTA only (IDH1 R132 in TARGET leukemia is not an HGG), and not when the
# methylation call is a control / low-purity tissue, whose lesions are unreliable

control_match_ids <- calls %>%
  filter(methyl_confidence == "high", methyl_control) %>%
  inner_join(select(hist, Kids_First_Biospecimen_ID, match_id),
             by = "Kids_First_Biospecimen_ID") %>%
  pull(match_id) %>%
  unique()
pbta_match_ids <- unique(hist$match_id[hist$cohort == "PBTA"])

fallback <- lesion_hits %>%
  filter(rule_id %in% rules$rule_id[rules$use_as_fallback],
         match_id %in% pbta_match_ids,
         !match_id %in% validated$match_id,
         !match_id %in% control_match_ids) %>%
  left_join(select(rules, rule_id, rule_order, fallback_requires_dx),
            by = "rule_id") %>%
  group_by(match_id) %>%
  arrange(rule_order, evidence, .by_group = TRUE) %>%
  summarise(OPC_molecular_subtype = first(OPC_molecular_subtype),
            lesion_evidence = paste(unique(evidence), collapse = "; "),
            requires_dx = first(fallback_requires_dx),
            .groups = "drop") %>%
  # a match_id can carry more than one diagnosis (e.g. cell-line rows), which
  # are joined with ";" as in the histologies
  left_join(hist %>%
              filter(sample_type == "Tumor") %>%
              group_by(match_id) %>%
              summarise(Kids_First_Participant_ID = first(Kids_First_Participant_ID),
                        sample_id = first(sample_id),
                        pathology_diagnosis = paste(sort(unique(na.omit(pathology_diagnosis))),
                                                    collapse = ";"),
                        .groups = "drop") %>%
              mutate(pathology_diagnosis = na_if(pathology_diagnosis, "")),
            by = "match_id") %>%
  mutate(dx_group = dx_to_groups(pathology_diagnosis, dx_groups))

# IDH-mutant gliomas: with 1p/19q codeletion the tumor is an oligodendroglioma
# (Oligo, IDH). Otherwise the grade comes from the diagnosis group, and a tumor
# with an oligodendroglioma diagnosis but no codeletion is held unclassified.
codeleted <- unique(lesion_hits$match_id[lesion_hits$rule_id == "OLIGO_1P19Q"])
arm_assessed <- unique(dna_ids$match_id[dna_ids$id %in%
                                          arm_values$Kids_First_Biospecimen_ID])

fallback <- fallback %>%
  mutate(
    is_idh = OPC_molecular_subtype == "HGG, IDH",
    is_codeleted = match_id %in% codeleted,
    oligo_dx = str_detect(replace_na(pathology_diagnosis, ""),
                          regex("oligodendroglioma", ignore_case = TRUE)),
    hold_reason = case_when(
      is_idh & !is_codeleted & oligo_dx & match_id %in% arm_assessed ~
        "oligodendroglioma diagnosis, 1p/19q not codeleted",
      is_idh & !is_codeleted & oligo_dx ~
        "oligodendroglioma diagnosis, 1p/19q not assessable (no WGS arm data)"
    ),
    OPC_molecular_subtype = case_when(
      # subtyped from the H3 variant itself
      OPC_molecular_subtype == "DMG, H3 K28" ~ "DMG, H3 K28-mutant",
      OPC_molecular_subtype == "DHG, H3 G35" ~ "DHG, H3 G35-mutant",
      is_idh & is_codeleted ~ "Oligo, IDH",
      is_idh & !oligo_dx & dx_group == "LGG" ~ "LGG, IDH",
      TRUE ~ OPC_molecular_subtype
    ),
    lesion_evidence = if_else(is_idh & is_codeleted,
                              paste0(lesion_evidence, "; 1p+19q codeleted"),
                              lesion_evidence),
    dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                    family_compat),
    evidence_level = "lesion only",
    # a call the diagnosis contradicts, or a held call, is not assigned; a
    # lesion that is not specific to the entity (BRAF V600E, CTNNB1 exon 3,
    # SMARCB1 loss) is only assigned when the diagnosis supports it
    assign = dx_concordance != "discordant" & is.na(hold_reason) &
      (!requires_dx | dx_concordance %in% c("concordant", "compatible"))
  ) %>%
  select(-is_idh, -is_codeleted, -oligo_dx, -requires_dx) %>%
  left_join(tp53, by = "match_id") %>%
  annotate_genes(annotation_genes, snv, dna_ids, cnv_segments) %>%
  annotate_expression(expression_annotations, expression_long, rna_ids) %>%
  add_h3_status(h3_hits, dna_ids) %>%
  relabel_egfr_only() %>%
  # fallback DMG tumors are subtyped from the H3 K27M/I variant itself
  add_h3k27_mechanism(h3_mutated = TRUE,
                      has_dna = .$match_id %in% dna_ids$match_id,
                      has_rna = .$match_id %in% rna_ids$match_id)

write_tsv(fallback, file.path(results_dir, "lesion_fallback_subtypes.tsv"))

### H3 variants in tumors that are not DMG/DHG --------------------------------
# Includes tumors with no subtype, e.g. where the methylation call is a control
# (low purity) and the lesion fallback is therefore not applied.
methyl_by_match <- calls %>%
  filter(methyl_confidence == "high") %>%
  inner_join(select(hist, Kids_First_Biospecimen_ID, match_id),
             by = "Kids_First_Biospecimen_ID") %>%
  distinct(match_id, methyl_class, methyl_control)

h3_variant_review <- h3_hits %>%
  filter(match_id %in% pbta_match_ids) %>%
  group_by(match_id) %>%
  summarise(h3_variant = paste(sort(unique(evidence)), collapse = "; "),
            .groups = "drop") %>%
  left_join(distinct(validated, match_id, methylation_subtype = OPC_molecular_subtype),
            by = "match_id") %>%
  left_join(select(fallback, match_id, fallback_subtype = OPC_molecular_subtype,
                   fallback_assigned = assign),
            by = "match_id") %>%
  left_join(methyl_by_match, by = "match_id") %>%
  left_join(hist %>%
              filter(sample_type == "Tumor") %>%
              group_by(match_id) %>%
              summarise(sample_id = first(sample_id),
                        pathology_diagnosis = paste(sort(unique(na.omit(pathology_diagnosis))),
                                                    collapse = ";"),
                        .groups = "drop"),
            by = "match_id") %>%
  mutate(
    reason = case_when(
      !is.na(methylation_subtype) & !str_detect(methylation_subtype, "^(DMG|DHG)") ~
        "methylation subtype is not DMG/DHG",
      methyl_control %in% TRUE ~
        "methylation is a control / low-purity call, so the lesion fallback is not applied",
      !is.na(fallback_assigned) & !fallback_assigned ~
        "lesion call held by the diagnosis",
      is.na(methylation_subtype) & is.na(fallback_subtype) ~ "no subtype assigned"
    )
  ) %>%
  filter(!is.na(reason)) %>%
  arrange(reason, sample_id)

write_tsv(h3_variant_review, file.path(results_dir, "h3_variant_review.tsv"))

### DMG from the free-text diagnosis ------------------------------------------
# PBTA tumors whose free-text diagnosis names a diffuse midline glioma or an H3
# K27 alteration, with no methylation call and no lesion-based subtype. The
# diagnosis is the only evidence, so the label is the plain "DMG, H3 K28" (not
# "-mutant", which needs the variant). "DIPG" alone is a clinical term and is
# not used. Tumors whose methylation is a control are left out, as for the
# lesion fallback.
dmg_free_text <- "diffuse midline|\\bdmg\\b|h3.?k27|h3 k28"

tumor_text <- hist %>%
  filter(sample_type == "Tumor", cohort == "PBTA", !is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(Kids_First_Participant_ID = first(Kids_First_Participant_ID),
            sample_id = first(sample_id),
            pathology_diagnosis = paste(sort(unique(na.omit(pathology_diagnosis))),
                                        collapse = ";"),
            free_text = paste(sort(unique(na.omit(pathology_free_text_diagnosis))),
                              collapse = " | "),
            .groups = "drop") %>%
  mutate(pathology_diagnosis = na_if(pathology_diagnosis, ""),
         free_text = na_if(free_text, ""))

free_text_dmg <- tumor_text %>%
  filter(str_detect(str_to_lower(replace_na(free_text, "")), dmg_free_text),
         !match_id %in% c(validated$match_id, fallback$match_id[fallback$assign],
                          control_match_ids)) %>%
  mutate(OPC_molecular_subtype = "DMG, H3 K28",
         lesion_evidence = "free-text diagnosis",
         evidence_level = "free text",
         dx_group = dx_to_groups(pathology_diagnosis, dx_groups),
         dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                         family_compat)) %>%
  left_join(tp53, by = "match_id") %>%
  annotate_genes(annotation_genes, snv, dna_ids, cnv_segments) %>%
  annotate_expression(expression_annotations, expression_long, rna_ids) %>%
  add_h3_status(h3_hits, dna_ids) %>%
  relabel_egfr_only() %>%
  add_h3k27_mechanism(h3_mutated = FALSE,
                      has_dna = .$match_id %in% dna_ids$match_id,
                      has_rna = .$match_id %in% rna_ids$match_id)

write_tsv(free_text_dmg, file.path(results_dir, "free_text_dmg.tsv"))

### Lesion-negative HGG -------------------------------------------------------
# PBTA tumors with an HGG diagnosis, a DNA assay, and no methylation call or
# lesion rule hit: no H3 or IDH variant means "HGG, H3 wildtype" by exclusion,
# which keeps the legacy catch-all without claiming a methylation subtype.
h3_idh_rule_ids <- c(h3_hits$rule_id, "HGG_IDH1", "HGG_IDH2", "OLIGO_IDH1",
                     "OLIGO_IDH2", "DMG_H3K28", "DHG_H3G35", "DMG_EGFR_H3K28")
has_any_rule_hit <- unique(lesion_hits$match_id)

lesion_negative_hgg <- hist %>%
  filter(sample_type == "Tumor", cohort == "PBTA", !is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(Kids_First_Participant_ID = first(Kids_First_Participant_ID),
            sample_id = first(sample_id),
            pathology_diagnosis = paste(sort(unique(na.omit(pathology_diagnosis))),
                                        collapse = ";"),
            .groups = "drop") %>%
  mutate(pathology_diagnosis = na_if(pathology_diagnosis, ""),
         dx_group = dx_to_groups(pathology_diagnosis, dx_groups)) %>%
  filter(str_detect(replace_na(dx_group, ""), "HGG"),
         match_id %in% dna_ids$match_id,
         !match_id %in% c(validated$match_id, fallback$match_id,
                          control_match_ids, has_any_rule_hit,
                          free_text_dmg$match_id)) %>%
  mutate(OPC_molecular_subtype = "HGG, H3 wildtype",
         lesion_evidence = "no H3 or IDH variant",
         evidence_level = "lesion negative",
         dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                         family_compat)) %>%
  left_join(tp53, by = "match_id") %>%
  annotate_genes(annotation_genes, snv, dna_ids, cnv_segments) %>%
  add_h3_status(h3_hits, dna_ids)

write_tsv(lesion_negative_hgg, file.path(results_dir, "lesion_negative_hgg.tsv"))

### Summary --------------------------------------------------------------------

message("\n== Methylation-based subtypes: ", nrow(validated))
print(count(validated, lesion_status))
print(count(validated, dx_concordance))
message("Flagged for review: ", sum(validated$needs_review))
message("Methylation overridden by a definitive lesion: ",
        sum(validated$methylation_overridden))
message("Methylation calls set aside (weak and contradicted by a definitive lesion): ",
        nrow(read_tsv(file.path(results_dir, "methylation_downgraded.tsv"),
                      show_col_types = FALSE)))
message("Lesion contradicts the methylation subtype: ",
        sum(!is.na(validated$lesion_contradiction)))
print(validated %>% filter(!is.na(lesion_contradiction)) %>%
        count(OPC_molecular_subtype, lesion_contradiction, sort = TRUE) %>%
        head(12))

message("\n== H3 variant but not DMG/DHG: ", nrow(h3_variant_review))
print(count(h3_variant_review, reason))
message("\n== DMG from the free-text diagnosis: ", nrow(free_text_dmg))
message("\n== Lesion-negative HGG, H3 wildtype: ", nrow(lesion_negative_hgg))
message("\n== Lesion-only fallback: ", nrow(fallback), " tumors without methylation")
print(count(fallback, OPC_molecular_subtype, assign, hold_reason))

# gaps in the lookup tables
gaps <- bind_rows(validated, fallback) %>%
  filter(dx_concordance == "subtype_family_unmapped" |
           str_detect(replace_na(dx_group, ""), "UNMAPPED"))
if (nrow(gaps) > 0) {
  message("\nSubtypes or diagnoses missing from the dx lookup tables:")
  print(distinct(gaps, OPC_molecular_subtype, pathology_diagnosis))
}
