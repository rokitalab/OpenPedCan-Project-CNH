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
         match_id, cohort, sample_type, experimental_strategy, RNA_library,
         primary_site, pathology_diagnosis, pathology_free_text_diagnosis)

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

# annoFuse RNA-Seq fusions plus the DGD panel fusions
fusion <- bind_rows(
  read_tsv(file.path(root_dir, "data", "fusion-annoFuse.tsv.gz"),
           show_col_types = FALSE) %>% select(Sample, FusionName),
  read_tsv(file.path(root_dir, "data", "fusion-dgd.tsv.gz"),
           show_col_types = FALSE) %>% select(Sample, FusionName)
) %>%
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

# manually reviewed subtypes decide which tumors the DMG annotations apply to
manual_subtype <- read_tsv(file.path(input_dir, "manual_subtype_review.tsv"),
                           show_col_types = FALSE) %>%
  transmute(sample_id, annotation_subtype = molecular_subtype)
with_manual <- function(df) left_join(df, manual_subtype, by = "sample_id")
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

# DNA: WGS, WXS and DNA panels; RNA: RNA-Seq and RNA panels (a targeted
# sequencing biospecimen with an RNA library, e.g. DGD fusions), as in the LGG
# module
dna_ids <- hist %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Sequencing"),
         is.na(RNA_library), !is.na(match_id)) %>%
  distinct(match_id, id = Kids_First_Biospecimen_ID)
rna_ids <- hist %>%
  filter(experimental_strategy == "RNA-Seq" |
           (experimental_strategy == "Targeted Sequencing" & !is.na(RNA_library)),
         !is.na(match_id)) %>%
  distinct(match_id, id = Kids_First_Biospecimen_ID)

# IHG RTK fusions from the in-frame, kinase-retained fusion summary (priority
# NTRK, ALK, ROS1, MET, as in the HGG module)
ihg <- read_tsv(file.path(root_dir, "data", "fusion_summary_lgg_hgg_foi.tsv"),
                show_col_types = FALSE) %>%
  pivot_longer(-Kids_First_Biospecimen_ID, names_to = "fusion", values_to = "n") %>%
  filter(n > 0) %>%
  inner_join(rna_ids, by = c("Kids_First_Biospecimen_ID" = "id"),
             relationship = "many-to-many") %>%
  mutate(term = case_when(
    str_detect(fusion, "NTRK") ~ "NTRK",
    str_detect(fusion, "ALK") ~ "ALK",
    str_detect(fusion, "ROS1") ~ "ROS1",
    str_detect(fusion, "MET") ~ "MET"
  )) %>%
  filter(!is.na(term)) %>%
  mutate(priority = match(term, c("NTRK", "ALK", "ROS1", "MET"))) %>%
  group_by(match_id) %>%
  filter(priority == min(priority)) %>%
  summarise(ihg_term = first(term),
            ihg_evidence = paste0(sort(unique(fusion)), " fusion", collapse = "; "),
            .groups = "drop")

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
  with_manual() %>%
  annotate_genes(annotation_genes, snv, dna_ids, cnv_segments) %>%
  annotate_expression(expression_annotations, expression_long, rna_ids) %>%
  add_h3_status(h3_hits, dna_ids) %>%
  relabel_egfr_only() %>%
  add_h3k27_mechanism(h3_mutated = .$lesion_status == "confirmed",
                      has_dna = .$match_id %in% dna_ids$match_id,
                      has_rna = .$match_id %in% rna_ids$match_id) %>%
  relabel_ihg(ihg) %>%
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
# an H3 K27/K28 M or I variant is specific to DMG, and a control methylation
# call is uninformative (low tumor purity), so these tumors are still subtyped
control_match_ids <- setdiff(
  control_match_ids,
  lesion_hits$match_id[lesion_hits$rule_id == "DMG_H3K28"])
# PBTA, plus PPTC xenograft models (e.g. ependymoma PDX), which are only
# assigned when the diagnosis is concordant or compatible
pbta_match_ids <- unique(hist$match_id[hist$cohort %in% c("PBTA", "PPTC")])
pptc_only_match_ids <- setdiff(unique(hist$match_id[hist$cohort == "PPTC"]),
                               unique(hist$match_id[hist$cohort == "PBTA"]))

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
                        free_text = paste(unique(na.omit(pathology_free_text_diagnosis)),
                                          collapse = " | "),
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
    # BRAF V600E is not specific, so PXA needs a PXA diagnosis in the free text
    is_pxa = OPC_molecular_subtype == "PXA",
    pxa_dx = str_detect(replace_na(free_text, ""),
                        regex("xanthoastrocytoma|\\bpxa\\b", ignore_case = TRUE)),
    pxa_prefix = case_when(
      str_detect(replace_na(dx_group, ""), "HGG") ~ "HGG",
      !is.na(dx_group) & !dx_group %in% c("UNINFORMATIVE", "UNMAPPED") ~ "LGG",
      free_text_grade_group(free_text) %in% "HGG" ~ "HGG",
      TRUE ~ "LGG"
    ),
    hold_reason = case_when(
      is_pxa & !pxa_dx ~ "BRAF V600E without a PXA diagnosis",
      OPC_molecular_subtype == "Rare CNS, DICER1 SARC" &
        str_detect(replace_na(dx_group, ""), "PINEAL") ~
        "DICER1 hotspot in pineoblastoma (PB, miRNA processing), not a DICER1 sarcoma",
      is_idh & !is_codeleted & oligo_dx & match_id %in% arm_assessed ~
        "oligodendroglioma diagnosis, 1p/19q not codeleted",
      is_idh & !is_codeleted & oligo_dx ~
        "oligodendroglioma diagnosis, 1p/19q not assessable (no WGS arm data)"
    ),
    OPC_molecular_subtype = case_when(
      # subtyped from the H3 variant itself
      OPC_molecular_subtype == "DMG, H3 K28" ~ "DMG, H3 K28-mutant",
      OPC_molecular_subtype == "DHG, H3 G35" ~ "DHG, H3 G35-mutant",
      is_pxa ~ paste0(pxa_prefix, ", PXA BRAF V600E"),
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
      (!requires_dx | dx_concordance %in% c("concordant", "compatible")) &
      (!match_id %in% pptc_only_match_ids |
         dx_concordance %in% c("concordant", "compatible"))
  ) %>%
  select(-is_idh, -is_codeleted, -oligo_dx, -requires_dx, -is_pxa, -pxa_dx,
         -pxa_prefix, -free_text) %>%
  left_join(tp53, by = "match_id") %>%
  with_manual() %>%
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
      fallback_assigned %in% TRUE ~ NA_character_,
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
            primary_site = first(na.omit(primary_site)),
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
  with_manual() %>%
  annotate_genes(annotation_genes, snv, dna_ids, cnv_segments) %>%
  annotate_expression(expression_annotations, expression_long, rna_ids) %>%
  add_h3_status(h3_hits, dna_ids) %>%
  relabel_egfr_only() %>%
  add_h3k27_mechanism(h3_mutated = FALSE,
                      has_dna = .$match_id %in% dna_ids$match_id,
                      has_rna = .$match_id %in% rna_ids$match_id)

# infant-type hemispheric glioma named in the free text, with the RTK fusion
# term when there is one
free_text_ihg <- tumor_text %>%
  filter(str_detect(str_to_lower(replace_na(free_text, "")),
                    "infant.?type hemispheric glioma"),
         !match_id %in% c(validated$match_id, fallback$match_id[fallback$assign],
                          control_match_ids, free_text_dmg$match_id)) %>%
  left_join(ihg, by = "match_id") %>%
  mutate(OPC_molecular_subtype = if_else(
           !is.na(ihg_term), paste0("HGG, IHG, ", ihg_term, "-altered"), "HGG, IHG"),
         lesion_evidence = coalesce(ihg_evidence, "free-text diagnosis"),
         evidence_level = "free text",
         dx_group = dx_to_groups(pathology_diagnosis, dx_groups),
         dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                         family_compat)) %>%
  select(-ihg_term, -ihg_evidence) %>%
  left_join(tp53, by = "match_id")

# adamantinomatous craniopharyngioma named in the free text, for tumors that
# methylation and the CTNNB1 / BRAF rules left without a subtype (a missing or
# low-confidence or control array: the free text is explicit about the histology
# and a control array only says the tissue was low purity). A high-confidence
# methylation call of another class stands and is listed in
# free_text_conflicts.tsv.
cranio_text <- "adamantinomato"
free_text_cranio <- tumor_text %>%
  filter(str_detect(str_to_lower(replace_na(free_text, "")), cranio_text),
         !match_id %in% c(validated$match_id, fallback$match_id[fallback$assign],
                          free_text_dmg$match_id, free_text_ihg$match_id)) %>%
  mutate(OPC_molecular_subtype = "CRANIO, ADAM",
         lesion_evidence = "free-text diagnosis",
         evidence_level = "free text",
         dx_group = dx_to_groups(pathology_diagnosis, dx_groups),
         dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                         family_compat)) %>%
  left_join(tp53, by = "match_id")

free_text_conflicts <- tumor_text %>%
  filter(str_detect(str_to_lower(replace_na(free_text, "")), cranio_text),
         match_id %in% validated$match_id) %>%
  inner_join(distinct(validated, match_id, methylation_subtype = OPC_molecular_subtype),
             by = "match_id") %>%
  filter(!str_detect(methylation_subtype, "^CRANIO")) %>%
  transmute(match_id, sample_id, free_text, methylation_subtype,
            free_text_conflict = "free-text diagnosis names adamantinomatous craniopharyngioma")
write_tsv(free_text_conflicts, file.path(results_dir, "free_text_conflicts.tsv"))

# neurocytoma (ported from molecular-subtyping-neurocytoma): the diagnosis
# decides, and the primary site separates central (ventricles, CNC) from
# extraventricular (EVN); a methylation call (CN, EVNCYT) takes precedence
neurocytoma_dx <- unlist(
  jsonlite::fromJSON(file.path(module_dir, "..", "molecular-subtyping-neurocytoma",
                               "neurocytoma-subset",
                               "neurocytoma_subtyping_path_dx_strings.json")))
free_text_neuro <- tumor_text %>%
  filter(pathology_diagnosis %in% neurocytoma_dx, !is.na(primary_site),
         !match_id %in% c(validated$match_id, fallback$match_id[fallback$assign],
                          control_match_ids, free_text_dmg$match_id,
                          free_text_ihg$match_id, free_text_cranio$match_id)) %>%
  mutate(OPC_molecular_subtype = if_else(primary_site == "Ventricles", "CNC", "EVN"),
         lesion_evidence = paste0("neurocytoma diagnosis, primary site ", primary_site),
         evidence_level = "diagnosis and primary site",
         dx_group = dx_to_groups(pathology_diagnosis, dx_groups),
         dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                         family_compat)) %>%
  left_join(tp53, by = "match_id")

free_text_dmg <- bind_rows(free_text_dmg, free_text_ihg, free_text_cranio,
                           free_text_neuro) %>%
  select(-any_of("primary_site"))

write_tsv(free_text_dmg, file.path(results_dir, "free_text_subtypes.tsv"))

### Embryonal: ETMR and CNS embryonal, NOS ------------------------------------
# Ported from molecular-subtyping-embryonal. Cohort: the embryonal diagnoses
# (and CNS neuroblastoma strings) of that module. ETMR, C19MC-altered is a
# TTYH1 fusion (5' partner) with LIN28A overexpression or a C19MC amplification;
# ETMR, NOS is an ETMR diagnosis without that, or LIN28A overexpression without
# TTYH1; the rest of the cohort is CNS Embryonal, NOS. LIN28A is the z-score
# within RNA library type (the module used z > 3 within the embryonal subset).
emb_strings <- jsonlite::fromJSON(file.path(input_dir, "embryonal_path_dx_strings.json"))
etmr_text <- "etmr|embryonal tumor with abundant neuropil and true rosettes|embryonal tumor with multilayer rosettes"

# a focal chr19 amplification overlapping the C19MC miRNA cluster (hg38)
c19mc_ids <- read_tsv(file.path(scratch_dir, "chr19_gain_segments.tsv"),
                      show_col_types = FALSE) %>%
  filter(loc.start <= 53767234, loc.end >= 53664934, seg.mean > 1,
         loc.end - loc.start < 2e6) %>%
  distinct(id = Kids_First_Biospecimen_ID) %>%
  inner_join(dna_ids, by = "id", relationship = "many-to-many") %>%
  distinct(match_id)

ttyh1 <- fusion %>%
  filter(str_detect(FusionName, "^TTYH1--")) %>%
  inner_join(rna_ids, by = c("Kids_First_Biospecimen_ID" = "id"),
             relationship = "many-to-many") %>%
  group_by(match_id) %>%
  summarise(ttyh1_fusions = paste(sort(unique(FusionName)), collapse = "; "),
            .groups = "drop")

lin28a <- expression_long %>%
  filter(gene == "LIN28A") %>%
  inner_join(rna_ids, by = c("Kids_First_Biospecimen_ID" = "id"),
             relationship = "many-to-many") %>%
  group_by(match_id) %>%
  summarise(lin28a_z = max(z, na.rm = TRUE), .groups = "drop")

embryonal <- hist %>%
  filter(sample_type == "Tumor", cohort %in% c("PBTA", "PPTC"), !is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(Kids_First_Participant_ID = first(Kids_First_Participant_ID),
            sample_id = first(sample_id),
            primary_site = paste(sort(unique(na.omit(primary_site))), collapse = ";"),
            pathology_diagnosis = paste(sort(unique(na.omit(pathology_diagnosis))), collapse = ";"),
            free_text = paste(sort(unique(na.omit(pathology_free_text_diagnosis))), collapse = " | "),
            .groups = "drop") %>%
  mutate(
    dx_l = str_to_lower(pathology_diagnosis),
    txt = str_to_lower(free_text),
    in_cohort =
      (str_detect(pathology_diagnosis,
                  paste(str_replace_all(emb_strings$include_path_dx, "([()])", "\\\\\\1"),
                        collapse = "|")) |
         str_detect(txt, paste(str_replace_all(emb_strings$include_free_text,
                                               "([()])", "\\\\\\1"), collapse = "|"))) &
        !str_detect(txt, paste(emb_strings$exclude_free_text, collapse = "|")) |
      # neuroblastoma diagnoses in the CNS
      (pathology_diagnosis == "Neuroblastoma" &
         !str_detect(primary_site, "Other locations NOS") &
         !str_detect(txt, "peripheral|metastatic"))
  ) %>%
  filter(in_cohort,
         !match_id %in% c(validated$match_id, fallback$match_id[fallback$assign],
                          control_match_ids, free_text_dmg$match_id)) %>%
  left_join(ttyh1, by = "match_id") %>%
  left_join(lin28a, by = "match_id") %>%
  mutate(
    has_ttyh1 = !is.na(ttyh1_fusions),
    lin28a_high = replace_na(lin28a_z >= 3, FALSE),
    c19mc_amp = match_id %in% c19mc_ids$match_id,
    etmr_dx = pathology_diagnosis == "Embryonal Tumor with Multilayered Rosettes (ETMR)" |
      str_detect(txt, etmr_text),
    OPC_molecular_subtype = case_when(
      (lin28a_high | c19mc_amp) & has_ttyh1 ~ "ETMR, C19MC-altered",
      etmr_dx | (lin28a_high & !has_ttyh1) ~ "ETMR, NOS",
      TRUE ~ "CNS Embryonal, NOS"
    ),
    lesion_evidence = pmap_chr(
      list(ttyh1_fusions, c19mc_amp, lin28a_high, lin28a_z),
      function(f, a, h, z) {
        e <- c(if (!is.na(f)) paste0(f, " fusion"),
               if (a) "C19MC amplification",
               if (h) paste0("LIN28A overexpression (z ", round(z, 1), ")"))
        if (length(e) == 0) NA_character_ else paste(e, collapse = "; ")
      }
    ),
    evidence_level = if_else(OPC_molecular_subtype == "CNS Embryonal, NOS",
                             "embryonal diagnosis", "fusion, copy number, expression"),
    dx_group = dx_to_groups(pathology_diagnosis, dx_groups),
    dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group, family_compat)
  ) %>%
  select(-dx_l, -txt, -in_cohort)

write_tsv(embryonal, file.path(results_dir, "embryonal_subtypes.tsv"))

### Lesion-negative HGG -------------------------------------------------------
# PBTA tumors with an HGG diagnosis, a DNA assay, and no methylation call or
# lesion rule hit: no H3 or IDH variant means "HGG, H3 wildtype" by exclusion,
# which keeps the legacy catch-all without claiming a methylation subtype.
h3_idh_rule_ids <- c(h3_hits$rule_id, "HGG_IDH1", "HGG_IDH2", "OLIGO_IDH1",
                     "OLIGO_IDH2", "DMG_H3K28", "DHG_H3G35", "DMG_EGFR_H3K28")
# only an H3 or IDH hit rules out the wildtype call; a hit in another rule
# (TP53, NF1, RB1, ...) is not specific and does not subtype the tumor
h3_idh_hit_ids <- unique(lesion_hits$match_id[lesion_hits$rule_id %in% h3_idh_rule_ids])

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
                          control_match_ids, h3_idh_hit_ids,
                          free_text_dmg$match_id)) %>%
  mutate(OPC_molecular_subtype = "HGG, H3 wildtype",
         lesion_evidence = "no H3 or IDH variant",
         evidence_level = "lesion negative",
         dx_concordance = dx_concordance(OPC_molecular_subtype, dx_group,
                                         family_compat)) %>%
  left_join(tp53, by = "match_id") %>%
  with_manual() %>%
  annotate_genes(annotation_genes, snv, dna_ids, cnv_segments) %>%
  annotate_expression(expression_annotations, expression_long, rna_ids) %>%
  add_h3_status(h3_hits, dna_ids) %>%
  add_h3k27_mechanism(h3_mutated = FALSE,
                      has_dna = .$match_id %in% dna_ids$match_id,
                      has_rna = .$match_id %in% rna_ids$match_id)

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
message("\n== Subtypes from the free-text diagnosis (DMG, IHG): ", nrow(free_text_dmg))
message("\n== Embryonal (ETMR, CNS Embryonal NOS): ", nrow(embryonal))
print(count(embryonal, OPC_molecular_subtype))
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
