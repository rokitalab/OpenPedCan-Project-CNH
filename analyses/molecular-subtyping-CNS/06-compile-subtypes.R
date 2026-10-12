################################################################################
# 06-compile-subtypes.R
#
# One row per match_id with the module's subtype and its evidence.
#
# Subtype precedence:
#   1. high-confidence methylation (03-validate-subtypes.R)
#   2. lesion-only fallback (03-validate-subtypes.R), when assigned
#   3. LGG lesion-level call, for tumors with no methylation or fallback call
#
# `lesion_refinement` adds the LGG lesion-level subtype (KIAA1549-BRAF, BRAF
# V600E, NF1, ...) on top of the methylation family when LGG made a specific
# call. TP53 columns are annotations and are not part of the subtype.
#
# Usage: Rscript --vanilla 06-compile-subtypes.R
################################################################################

suppressPackageStartupMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS")
results_dir <- file.path(module_dir, "results")

# free-text evidence columns are often entirely empty, so force them to character
evidence_types <- cols(tp53_snv_evidence = col_character(),
                       tp53_cnv_evidence = col_character(),
                       egfr_variants = col_character(),
                       .default = col_guess())

# chordoma: SMARCB1 focal copy status and expression, as annotations
smarcb1 <- read_tsv(
  file.path(root_dir, "analyses", "molecular-subtyping-chordoma", "results",
            "chordoma_smarcb1_status.tsv"),
  col_types = cols(.default = col_character())
) %>%
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(
    smarcb1_sample_id = first(sample_id),
    smarcb1_participant = first(Kids_First_Participant_ID),
    smarcb1_focal_status = first(na.omit(focal_SMARCB1_status)),
    smarcb1_expression = first(na.omit(SMARCB1_expression)),
    .groups = "drop"
  )

# manual review (input/manual_subtype_review.tsv): a subtype set by a pathology
# report or clinical review. It has the highest precedence; the automated call is
# kept in automated_subtype and disagreements with it are flagged.
manual_review <- read_tsv(file.path(module_dir, "input", "manual_subtype_review.tsv"),
                          show_col_types = FALSE) %>%
  transmute(sample_id, manual_subtype = molecular_subtype,
            manual_source = source, manual_note = note)
stopifnot(!anyDuplicated(manual_review$sample_id))

# match_id for each reviewed sample (tumor rows of histologies-base)
tumor_ids <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                      guess_max = 100000, show_col_types = FALSE) %>%
  filter(sample_type == "Tumor", !is.na(match_id),
         sample_id %in% manual_review$sample_id) %>%
  distinct(sample_id, match_id, Kids_First_Participant_ID) %>%
  group_by(sample_id) %>%
  slice(1) %>%
  ungroup()

# umbrella groups for subtypes (input/subtype_groups.tsv)
subtype_groups <- read_tsv(file.path(module_dir, "input", "subtype_groups.tsv"),
                           show_col_types = FALSE)

umbrella_group <- function(subtype) {
  if (is.na(subtype)) return(NA_character_)
  hit <- subtype_groups %>% filter(str_detect(subtype, subtype_regex)) %>% slice(1)
  if (nrow(hit) == 0) NA_character_ else hit$subtype_group
}

# LGG lesion-level subtype from 04-lgg-lesions.R (ported from the LGAT module)
lgg <- read_tsv(file.path(results_dir, "lgg_lesions.tsv"),
                show_col_types = FALSE) %>%
  transmute(match_id, lgg_subtype = lgg_lesion_label, lgg_alterations,
            lgg_cohort, lgg_specific = lgg_lesion_status == "specific",
            lgg_status = lgg_lesion_status, lgg_has_dna = has_dna,
            lgg_has_rna = has_rna)

# medulloblastoma: expression classifier, concordance with methylation, and
# the SHH subtype (05-validate-mb.R)
mb <- read_tsv(file.path(results_dir, "mb_validation.tsv"),
               guess_max = 10000, show_col_types = FALSE) %>%
  filter(!is.na(match_id)) %>%
  select(match_id, mb_sample_id = sample_id, mb_expression_group,
         mb_concordance, mb_shh_subtype, mb_shh_source)

# one methylation call per tumor: highest score wins when several biospecimens
methyl <- read_tsv(file.path(results_dir, "methylation_subtype_validation.tsv"),
                   guess_max = 10000, col_types = evidence_types) %>%
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  arrange(desc(methyl_score), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(subtype_source = if_else(methylation_overridden %in% TRUE,
                                  "lesion overrides methylation", "methylation"))

fallback <- read_tsv(file.path(results_dir, "lesion_fallback_subtypes.tsv"),
                     guess_max = 10000, col_types = evidence_types) %>%
  filter(assign, !match_id %in% methyl$match_id) %>%
  mutate(subtype_source = "lesion fallback")

free_text <- read_tsv(file.path(results_dir, "free_text_subtypes.tsv"),
                      guess_max = 10000, col_types = evidence_types) %>%
  filter(!match_id %in% c(methyl$match_id, fallback$match_id)) %>%
  mutate(subtype_source = if_else(evidence_level == "diagnosis and primary site",
                                  "pathology diagnosis", "pathology free text"))

# free-text diagnoses that disagree with a methylation-based subtype
free_text_conflicts <- read_tsv(file.path(results_dir, "free_text_conflicts.tsv"),
                                show_col_types = FALSE) %>%
  distinct(match_id, free_text_conflict)

embryonal <- read_tsv(file.path(results_dir, "embryonal_subtypes.tsv"),
                      guess_max = 10000, col_types = evidence_types) %>%
  filter(!match_id %in% c(methyl$match_id, fallback$match_id)) %>%
  mutate(subtype_source = "embryonal rules")

lesion_negative <- read_tsv(file.path(results_dir, "lesion_negative_hgg.tsv"),
                            guess_max = 10000, col_types = evidence_types) %>%
  filter(!match_id %in% c(methyl$match_id, fallback$match_id)) %>%
  mutate(subtype_source = "lesion negative")

shared <- c("match_id", "Kids_First_Participant_ID", "sample_id",
            "OPC_molecular_subtype", "subtype_source", "evidence_level",
            "lesion_evidence", "pathology_diagnosis", "dx_group",
            "dx_concordance")

compiled <- bind_rows(
  methyl %>% select(all_of(shared), methyl_class, methyl_score, lesion_status,
                    lesion_contradiction,
                    methyl_conflict, starts_with("tp53_"),
                    ends_with("_mutation_status"), ends_with("_variants"),
                    ends_with("_amplification"), ends_with("_activating_alteration"),
                    ends_with("_tpm"), ends_with("_z"), ends_with("_expression_status"),
                    starts_with("h3_"), starts_with("h3k27_")),
  embryonal %>% select(all_of(shared)),
  free_text %>% select(all_of(shared), starts_with("tp53_"),
                       ends_with("_mutation_status"), ends_with("_variants"),
                       ends_with("_amplification"),
                       ends_with("_activating_alteration"),
                       ends_with("_tpm"), ends_with("_z"),
                       ends_with("_expression_status"),
                       starts_with("h3_"), starts_with("h3k27_")),
  lesion_negative %>% select(all_of(shared), starts_with("tp53_"),
                             ends_with("_mutation_status"), ends_with("_variants"),
                             ends_with("_amplification"),
                             ends_with("_activating_alteration"),
                             ends_with("_tpm"), ends_with("_z"),
                             ends_with("_expression_status"),
                             starts_with("h3_"), starts_with("h3k27_")),
  fallback %>% select(all_of(shared), starts_with("tp53_"),
                      ends_with("_mutation_status"), ends_with("_variants"),
                    ends_with("_amplification"), ends_with("_activating_alteration"),
                    ends_with("_tpm"), ends_with("_z"), ends_with("_expression_status"),
                    starts_with("h3_"), starts_with("h3k27_"))
) %>%
  # LGG lesions annotate every tumor already here; tumors only the LGG cohort
  # brings in are added as rows (the lesions are computed for all PBTA tumors,
  # so they are not all added)
  left_join(lgg, by = "match_id") %>%
  (function(d) bind_rows(d, filter(lgg, lgg_cohort, !match_id %in% d$match_id))) %>%
  full_join(mb, by = "match_id") %>%
  # full join: chordoma tumors with no other source can still be subtyped
  # from focal SMARCB1 loss
  full_join(smarcb1, by = "match_id") %>%
  left_join(free_text_conflicts, by = "match_id") %>%
  mutate(
    sample_id = coalesce(sample_id, mb_sample_id, smarcb1_sample_id),
    Kids_First_Participant_ID = coalesce(Kids_First_Participant_ID,
                                         smarcb1_participant),
    # the LGG lesion call is the subtype only for the LGG cohort, and only when
    # nothing else assigned one
    # expression classifier group, when it names a single group
    mb_expression_single = if_else(str_detect(replace_na(mb_expression_group, ""), "\\|"),
                                   NA_character_, mb_expression_group),
    subtype_source = case_when(
      !is.na(subtype_source) ~ subtype_source,
      lgg_specific %in% TRUE & lgg_cohort %in% TRUE ~ "LGG lesion",
      !is.na(mb_expression_single) ~ "expression classifier",
      # chordoma diagnosis (the chordoma module's cohort) with focal SMARCB1 loss
      smarcb1_focal_status %in% "loss" ~ "SMARCB1 focal loss",
      TRUE ~ NA_character_
    ),
    molecular_subtype_cns = case_when(
      !is.na(OPC_molecular_subtype) ~ OPC_molecular_subtype,
      subtype_source == "LGG lesion" ~ lgg_subtype,
      subtype_source == "expression classifier" ~ mb_expression_single,
      subtype_source == "SMARCB1 focal loss" ~ "CHDM, poorly differentiated"
    ),
    evidence_level = case_when(
      subtype_source == "pathology free text" ~ "free text",
      subtype_source %in% c("LGG lesion", "SMARCB1 focal loss") ~ "lesion only",
      subtype_source == "embryonal rules" ~ evidence_level,
      subtype_source == "expression classifier" ~ "expression only",
      TRUE ~ evidence_level
    ),
    methyl_conflict = replace_na(methyl_conflict, FALSE)
  ) %>%
  # reviewed tumors with no automated source still need a row
  (function(d) {
    missing <- manual_review %>%
      filter(!sample_id %in% d$sample_id) %>%
      inner_join(tumor_ids, by = "sample_id")
    bind_rows(d, missing %>% select(match_id, Kids_First_Participant_ID, sample_id))
  }) %>%
  # manual review overrides the automated subtype
  left_join(manual_review, by = "sample_id") %>%
  mutate(
    automated_subtype = molecular_subtype_cns,
    manual_conflict = !is.na(manual_subtype) & !is.na(automated_subtype) &
      manual_subtype != automated_subtype,
    subtype_source = if_else(!is.na(manual_subtype), "manual review", subtype_source),
    evidence_level = if_else(!is.na(manual_subtype), "manual review", evidence_level),
    molecular_subtype_cns = coalesce(manual_subtype, molecular_subtype_cns),
    # EZHIP overexpression is added to reviewed DMG labels as well
    molecular_subtype_cns = case_when(
      !str_detect(replace_na(molecular_subtype_cns, ""), "^DMG") ~ molecular_subtype_cns,
      !(ezhip_expression_status %in% "overexpressed") ~ molecular_subtype_cns,
      str_detect(molecular_subtype_cns, "EZHIP") ~ molecular_subtype_cns,
      molecular_subtype_cns == "DMG, H3 K28" ~ "DMG, H3 K28-altered, EZHIP",
      TRUE ~ paste0(molecular_subtype_cns, ", EZHIP")
    )
  ) %>%
  # lesion-level refinement from the LGG lesions, for the LGG / GNT / GNG / SEGA
  # families (the LGG cohort's own call is the subtype, shown here too)
  mutate(lesion_refinement = if_else(
    lgg_specific %in% TRUE &
      (str_detect(replace_na(molecular_subtype_cns, ""), "^(LGG|GNT|GNG|SEGA)") |
         subtype_source %in% "LGG lesion"),
    lgg_subtype, NA_character_
  )) %>%
  # what the lesion refinement could see: specific (a lesion term), none
  # detected (DNA and RNA available), or incomplete (a missing assay, so the
  # lesion is still to be determined, not absent)
  mutate(
    in_lgg_families = str_detect(replace_na(molecular_subtype_cns, ""),
                                 "^(LGG|GNT|GNG|SEGA)") |
      subtype_source %in% "LGG lesion",
    lesion_refinement_status = case_when(
      !in_lgg_families | is.na(lgg_status) ~ NA_character_,
      !is.na(lesion_refinement) ~ "specific",
      lgg_status == "none detected" ~ "none detected (DNA and RNA)",
      lgg_has_dna %in% FALSE & lgg_has_rna %in% FALSE ~ "incomplete: no DNA or RNA",
      lgg_has_dna %in% FALSE ~ "incomplete: no DNA",
      lgg_has_rna %in% FALSE ~ "incomplete: no RNA",
      TRUE ~ NA_character_
    )
  ) %>%
  # the supratentorial pilocytic astrocytoma classes carry their lesion term in
  # the subtype (e.g. GNG, PA ST, BRAF V600E); without a lesion the label stays
  # at the class and lesion_refinement_status says why
  mutate(molecular_subtype_cns = if_else(
    str_detect(replace_na(molecular_subtype_cns, ""), ", PA ST$") &
      !is.na(lesion_refinement),
    paste0(molecular_subtype_cns, ", ", str_remove(lesion_refinement, "^[^,]*, ")),
    molecular_subtype_cns
  )) %>%
  # umbrella group (e.g. the pediatric H3-wildtype HGG classes), unless an H3
  # variant was found, which contradicts the wildtype grouping
  mutate(subtype_group = map_chr(molecular_subtype_cns, umbrella_group)) %>%
  mutate(
    h3_group_conflict = !is.na(subtype_group) & h3_status %in% "mutant",
    subtype_group = case_when(
      h3_group_conflict ~ molecular_subtype_cns,
      !is.na(subtype_group) ~ subtype_group,
      TRUE ~ molecular_subtype_cns
    ),
    needs_review = is.na(manual_subtype) &
      (dx_concordance %in% "discordant" | methyl_conflict |
         h3_group_conflict | h3_variant_conflict %in% TRUE |
         mb_concordance %in% "discordant" | !is.na(lesion_contradiction) |
         !is.na(free_text_conflict))
  ) %>%
  # the genomic and expression alterations that confirm the subtype: the rule
  # evidence that matched (lesion_evidence), the H3 variant for -mutant labels,
  # the EGFR alteration for EGFR labels, EZHIP overexpression for EZHIP labels, and
  # the expression classifier for medulloblastoma
  mutate(
    # EGFR evidence for EGFR labels, per tumor
    egfr_variant_text = str_replace_all(egfr_activating_variants,
                                        "(p\\.[^;]+)", "EGFR \\1"),
    egfr_evidence = case_when(
      !str_detect(replace_na(molecular_subtype_cns, ""), "EGFR") ~ NA_character_,
      !is.na(egfr_variant_text) & egfr_amplification %in% "amplified" ~
        paste(egfr_variant_text, "EGFR amplification", sep = "; "),
      !is.na(egfr_variant_text) ~ egfr_variant_text,
      egfr_amplification %in% "amplified" ~ "EGFR amplification",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    subtype_alterations_confirmed = pmap_chr(
      list(
        if_else(lesion_status %in% "confirmed" |
                  subtype_source %in% c("lesion fallback", "lesion overrides methylation"),
                lesion_evidence, NA_character_),
        if_else(h3_status %in% "mutant" & str_detect(replace_na(molecular_subtype_cns, ""),
                                                     "-mutant"), h3_variant, NA_character_),
        egfr_evidence,
        if_else(str_detect(replace_na(molecular_subtype_cns, ""), "EZHIP"),
                "EZHIP overexpression", NA_character_),
        # expression evidence for medulloblastoma: the expression classifier
        # agrees with the methylation group
        if_else(mb_concordance %in% "concordant",
                paste0("expression classifier agrees (", mb_expression_group, ")"),
                NA_character_),
        if_else(subtype_source %in% "SMARCB1 focal loss", "SMARCB1 focal loss",
                NA_character_),
        # lesion-level evidence for the LGG / GNT / GNG / SEGA families
        if_else(!is.na(lesion_refinement), lgg_alterations, NA_character_)
      ),
      function(...) {
        x <- unlist(strsplit(na.omit(c(...)), "; ", fixed = TRUE))
        x <- unique(x[x != ""])
        if (length(x) == 0) NA_character_ else paste(sort(x), collapse = "; ")
      }
    ),
    alteration_confirmed = !is.na(subtype_alterations_confirmed)
  ) %>%
  # provenance for the release file (a fresh vocabulary, not the legacy Notes
  # strings, since the subtypes are being redone); subtype_review_flag is a
  # short coded reason
  mutate(
    subtype_review_flag = pmap_chr(
      list(is.na(manual_subtype) & dx_concordance %in% "discordant",
           is.na(manual_subtype) & methyl_conflict,
           is.na(manual_subtype) & !is.na(lesion_contradiction),
           is.na(manual_subtype) & h3_variant_conflict %in% TRUE,
           is.na(manual_subtype) & h3_group_conflict %in% TRUE,
           is.na(manual_subtype) & mb_concordance %in% "discordant",
           is.na(manual_subtype) & !is.na(free_text_conflict)),
      function(...) {
        flags <- c("diagnosis discordant", "classifier versions disagree",
                   "lesion contradicts methylation", "H3 variant outside DMG/DHG",
                   "H3 variant in wildtype group", "expression classifier discordant",
                   "free-text diagnosis differs from methylation")
        hit <- flags[c(...)]
        if (length(hit) == 0) NA_character_ else paste(hit, collapse = ";")
      }
    ),
    Notes = case_when(
      subtype_source == "manual review" ~
        "Subtype from manual review of pathology report",
      subtype_source == "methylation" ~ "Subtype from methylation",
      subtype_source == "lesion overrides methylation" ~
        "Subtype from lesion, overriding methylation",
      subtype_source %in% c("lesion fallback", "LGG lesion", "SMARCB1 focal loss") ~
        "Subtype from lesion",
      subtype_source == "lesion negative" ~
        "Subtype by exclusion, no H3 or IDH variant",
      subtype_source == "expression classifier" ~
        "Subtype from expression classifier",
      subtype_source == "pathology free text" ~
        "Subtype from pathology free text diagnosis",
      subtype_source == "pathology diagnosis" ~
        "Subtype from pathology diagnosis and primary site",
      subtype_source == "embryonal rules" &
        str_detect(replace_na(molecular_subtype_cns, ""), "^ETMR") ~
        "Subtype from fusion, copy number and expression",
      subtype_source == "embryonal rules" ~
        "Subtype by exclusion, embryonal diagnosis",
      TRUE ~ NA_character_
    )
  ) %>%
  select(match_id, Kids_First_Participant_ID, sample_id,
         molecular_subtype_cns, subtype_group, lesion_refinement,
         lesion_refinement_status, subtype_source,
         subtype_alterations_confirmed, alteration_confirmed, subtype_review_flag, Notes,
         automated_subtype, manual_conflict, manual_source, manual_note,
         evidence_level, lesion_status, lesion_evidence, lesion_contradiction,
         methyl_class,
         mb_expression_group, mb_concordance, mb_shh_subtype,
         smarcb1_focal_status, smarcb1_expression,
         methyl_score, methyl_conflict, pathology_diagnosis, dx_group,
         dx_concordance, needs_review, free_text_conflict, starts_with("tp53_"),
         ends_with("_mutation_status"), ends_with("_variants"),
         ends_with("_amplification"), ends_with("_activating_alteration"),
         ends_with("_tpm"), ends_with("_z"), ends_with("_expression_status"),
         starts_with("h3_"), starts_with("h3k27_")) %>%
  arrange(match_id)

stopifnot(!anyDuplicated(compiled$match_id))

### Cell lines versus the tumor --------------------------------------
# A derived cell line is subtyped from its own methylation and lesions (its own
# match_id). The tumors are the solid-tissue tumors of the same sample
# (or, failing that, of the same participant). Concordance: same = identical
# subtype, same group = the same umbrella group (e.g. HGG, H3 wildtype for an
# HGG, RTK tumor), different = another subtype, and no subtype on either side.
composition <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                        guess_max = 100000, show_col_types = FALSE) %>%
  filter(!is.na(match_id)) %>%
  distinct(match_id, composition)
is_line <- composition$match_id[composition$composition %in% "Derived Cell Line"]

tumors_by_line <- compiled %>%
  filter(!match_id %in% is_line, !is.na(molecular_subtype_cns),
         !is.na(sample_id), !is.na(Kids_First_Participant_ID)) %>%
  select(Kids_First_Participant_ID, sample_id, tumor_subtype = molecular_subtype_cns,
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
    tumor_concordance = case_when(
      !is_cell_line ~ NA_character_,
      is.na(molecular_subtype_cns) | is.na(tumor_subtype) ~ "no subtype to compare",
      molecular_subtype_cns == tumor_subtype ~ "same",
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
  filter(!match_id %in% is_line, !is.na(subtype_alterations_confirmed),
         !is.na(sample_id)) %>%
  select(sample_id, Kids_First_Participant_ID, tumor_subtype = molecular_subtype_cns,
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
  filter(!is.na(tumor_concordance)) %>%
  left_join(line_alterations, by = "match_id") %>%
  mutate(res = pmap(
    list(sample_id, Kids_First_Participant_ID, molecular_subtype_cns,
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
    molecular_subtype_cns = if_else(adopted, adopt_subtype, molecular_subtype_cns),
    subtype_group = if_else(adopted, adopt_subtype, subtype_group),
    subtype_source = if_else(adopted, "lesion shared with tumor", subtype_source),
    evidence_level = if_else(adopted, "lesion shared with tumor", evidence_level),
    subtype_alterations_confirmed = if_else(adopted, adopt_alterations,
                                            subtype_alterations_confirmed),
    alteration_confirmed = alteration_confirmed | adopted,
    Notes = if_else(adopted, "Subtype from lesion shared with the tumor", Notes),
    tumor_concordance = case_when(
      is.na(tumor_concordance) ~ NA_character_,
      !is.na(molecular_subtype_cns) & molecular_subtype_cns == tumor_subtype ~ "same",
      # one label is the other plus an annotation term (e.g. ", EGFR", "-mutant")
      !is.na(molecular_subtype_cns) & !is.na(tumor_subtype) &
        (str_starts(tumor_subtype, fixed(molecular_subtype_cns)) |
           str_starts(molecular_subtype_cns, fixed(tumor_subtype))) ~ "same group",
      adopted & subtype_group == tumor_group ~ "same group",
      TRUE ~ tumor_concordance
    )
  ) %>%
  select(-adopt_subtype, -adopt_alterations, -adopted)
message("Cell lines taking the tumor subtype from a shared lesion: ",
        sum(compiled$subtype_source %in% "lesion shared with tumor"))

# a line that still differs from its tumor is flagged for review
compiled <- compiled %>%
  mutate(
    needs_review = needs_review | tumor_concordance %in% "different",
    subtype_review_flag = case_when(
      tumor_concordance %in% "different" & is.na(subtype_review_flag) ~
        "cell line differs from tumor",
      tumor_concordance %in% "different" ~
        paste(subtype_review_flag, "cell line differs from tumor", sep = ";"),
      TRUE ~ subtype_review_flag
    )
  )

write_tsv(
  compiled %>% filter(!is.na(tumor_concordance)) %>%
    select(match_id, sample_id, molecular_subtype_cns, subtype_source,
           tumor_subtype, tumor_concordance, tumor_alterations_shared),
  file.path(results_dir, "cell_line_tumor_concordance.tsv")
)
write_tsv(compiled, file.path(results_dir, "molecular_subtype_cns.tsv"))

message("Tumors with a subtype: ",
        sum(!is.na(compiled$molecular_subtype_cns)), " of ", nrow(compiled))
print(count(compiled, subtype_source))
message("With lesion refinement: ", sum(!is.na(compiled$lesion_refinement)))
message("\nNotes")
print(count(compiled, Notes, sort = TRUE))
message("\nReview flags")
print(count(compiled, subtype_review_flag, sort = TRUE))
