################################################################################
# 07-compile-subtypes.R
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
# Usage: Rscript --vanilla 07-compile-subtypes.R
################################################################################

suppressPackageStartupMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
results_dir <- file.path(module_dir, "results")
source(file.path(module_dir, "utils", "cell-lines.R"))
source(file.path(module_dir, "utils", "provenance.R"))

# free-text evidence columns are often entirely empty, so force them to character
evidence_types <- cols(tp53_snv_evidence = col_character(),
                       tp53_cnv_evidence = col_character(),
                       egfr_variants = col_character(),
                       .default = col_guess())

# chordoma: SMARCB1 focal copy status and expression, as annotations. Focal
# SMARCB1 loss is the lesion of poorly differentiated chordoma. The status is
# "loss" for the WGS/WXS biospecimens with a focal SMARCB1 loss and "neutral" for
# the other WGS/WXS biospecimens; the expression is the SMARCB1 TPM of the RNA-seq
# biospecimens; both are collapsed per tumor (match_id).
source(file.path(module_dir, "utils", "duckdb-connect.R"))
source(file.path(module_dir, "utils", "release-tables.R"))
hist_base <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                      guess_max = 100000, show_col_types = FALSE)
# tumors whose methylation is a likely or possible mismatch of the specimen
# (results/methylation_mismatch_calls.tsv, from 03-validate-subtypes.R: found from the
# lesions and the diagnosis, and the reviewer decisions of
# input/possible_methylation_mixups.tsv)
mismatch <- read_tsv(file.path(results_dir, "methylation_mismatch_calls.tsv"),
                     show_col_types = FALSE, col_types = cols(.default = "c"))
stopifnot(all(mismatch$mismatch_call %in% c("likely", "possible", "kept")),
          !anyDuplicated(mismatch$sample_id))

# diagnosis discordances reviewed with the pathology report (or an obsolete
# diagnosis): the methylation call is accepted and the flag is cleared
reviewed <- read_tsv(file.path(module_dir, "input", "reviewed_diagnosis_discordance.tsv"),
                     show_col_types = FALSE, col_types = cols(.default = "c"))
stopifnot(!anyDuplicated(reviewed$sample_id))

chordoma_dx <- jsonlite::fromJSON(
  file.path(module_dir, "input", "path_dx_strings.json"))$chordoma$exact_path_dx
chordoma <- hist_base %>%
  filter(pathology_diagnosis == chordoma_dx) %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, match_id,
         sample_id, experimental_strategy)

smarcb1_con <- connect_duckdb()
smarcb1_loss <- DBI::dbGetQuery(
  smarcb1_con,
  sprintf("SELECT DISTINCT biospecimen_id FROM %s
           WHERE gene_symbol = 'SMARCB1' AND status = 'loss'",
          release_source("consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only",
                         "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz",
                         root_dir)))$biospecimen_id
DBI::dbDisconnect(smarcb1_con, shutdown = TRUE)

expression_all <- readRDS(file.path(root_dir, "data",
                                    "gene-expression-rsem-tpm-collapsed.rds"))
smarcb1_tpm <- unlist(expression_all["SMARCB1", ])
rm(expression_all)

smarcb1_status <- chordoma %>%
  mutate(
    focal_SMARCB1_status = case_when(
      !experimental_strategy %in% c("WGS", "WXS") ~ NA_character_,
      Kids_First_Biospecimen_ID %in% smarcb1_loss ~ "loss",
      TRUE ~ "neutral"),
    SMARCB1_expression = unname(smarcb1_tpm[Kids_First_Biospecimen_ID]))

smarcb1_by_tumor <- smarcb1_status %>%
  group_by(match_id) %>%
  summarise(
    focal_SMARCB1_status = paste(unique(na.omit(focal_SMARCB1_status)), collapse = ", "),
    SMARCB1_expression = paste(unique(na.omit(SMARCB1_expression)), collapse = ", "),
    .groups = "drop")

smarcb1 <- smarcb1_status %>%
  distinct(match_id, .keep_all = TRUE) %>%
  select(match_id, smarcb1_sample_id = sample_id,
         smarcb1_participant = Kids_First_Participant_ID) %>%
  left_join(smarcb1_by_tumor, by = "match_id") %>%
  transmute(match_id, smarcb1_sample_id, smarcb1_participant,
            smarcb1_focal_status = na_if(focal_SMARCB1_status, ""),
            smarcb1_expression = na_if(SMARCB1_expression, ""))
write_tsv(smarcb1, file.path(results_dir, "chordoma_smarcb1_status.tsv"))

# manual review (input/manual_subtype_review.tsv): a subtype set by a pathology
# report or clinical review. It has the highest precedence; the automated call is
# kept in automated_subtype and disagreements with it are flagged.
manual_review <- read_manual_review(module_dir)

# match_id for each reviewed sample (tumor rows of histologies-base)
tumor_ids <- reviewed_tumor_ids(manual_review, root_dir)

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

# solid tumors (06-subtype-solid-tumors.R): neuroblastic tumors by MYCN status.
# The subtype from MYCN replaces a methylation-class label for these tumors.
solid <- read_tsv(file.path(results_dir, "solid_subtypes.tsv"),
                  show_col_types = FALSE)
solid_called <- solid %>% filter(!is.na(molecular_subtype))

shared <- c("match_id", "Kids_First_Participant_ID", "sample_id",
            "OPC_molecular_subtype", "subtype_source", "evidence_level",
            "lesion_evidence", "pathology_diagnosis", "dx_group",
            "dx_concordance")

compiled <- bind_rows(
  solid_called %>%
    transmute(match_id, Kids_First_Participant_ID, sample_id,
              OPC_molecular_subtype = molecular_subtype, subtype_source,
              evidence_level, lesion_evidence = subtype_alterations_confirmed,
              pathology_diagnosis, dx_group, dx_concordance),
  methyl %>% filter(!match_id %in% solid_called$match_id) %>%
    select(all_of(shared), methyl_class, methyl_score, lesion_status,
                    lesion_contradiction,
                    methyl_conflict, starts_with("tp53_"),
                    ends_with("_mutation_status"), ends_with("_variants"),
                    ends_with("_amplification"), ends_with("_activating_alteration"),
                    ends_with("_tpm"), ends_with("_z"), ends_with("_expression_status"),
                    starts_with("h3_"), starts_with("h3k27_")),
  embryonal %>% filter(!match_id %in% solid_called$match_id) %>%
    select(all_of(shared)),
  free_text %>% filter(!match_id %in% solid_called$match_id) %>%
    select(all_of(shared), starts_with("tp53_"),
                       ends_with("_mutation_status"), ends_with("_variants"),
                       ends_with("_amplification"),
                       ends_with("_activating_alteration"),
                       ends_with("_tpm"), ends_with("_z"),
                       ends_with("_expression_status"),
                       starts_with("h3_"), starts_with("h3k27_")),
  lesion_negative %>% filter(!match_id %in% solid_called$match_id) %>%
    select(all_of(shared), starts_with("tp53_"),
                             ends_with("_mutation_status"), ends_with("_variants"),
                             ends_with("_amplification"),
                             ends_with("_activating_alteration"),
                             ends_with("_tpm"), ends_with("_z"),
                             ends_with("_expression_status"),
                             starts_with("h3_"), starts_with("h3k27_")),
  fallback %>% filter(!match_id %in% solid_called$match_id) %>%
    select(all_of(shared), starts_with("tp53_"),
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
  left_join(solid %>% select(match_id, solid_review_flag = subtype_review_flag,
                             mycn_cn_status, mycn_tpm, mycn_text_status),
            by = "match_id") %>%
  (function(d) bind_rows(d, solid %>%
     filter(!is.na(subtype_review_flag), !match_id %in% d$match_id) %>%
     transmute(match_id, Kids_First_Participant_ID, sample_id,
               solid_review_flag = subtype_review_flag, mycn_cn_status,
               mycn_tpm, mycn_text_status, pathology_diagnosis))) %>%
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
    molecular_subtype = case_when(
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
    automated_subtype = molecular_subtype,
    manual_conflict = !is.na(manual_subtype) & !is.na(automated_subtype) &
      manual_subtype != automated_subtype,
    subtype_source = if_else(!is.na(manual_subtype), "manual review", subtype_source),
    evidence_level = if_else(!is.na(manual_subtype), "manual review", evidence_level),
    molecular_subtype = coalesce(manual_subtype, molecular_subtype),
    # EZHIP overexpression is added to reviewed DMG labels as well
    molecular_subtype = case_when(
      !str_detect(replace_na(molecular_subtype, ""), "^DMG") ~ molecular_subtype,
      !(ezhip_expression_status %in% "overexpressed") ~ molecular_subtype,
      str_detect(molecular_subtype, "EZHIP") ~ molecular_subtype,
      molecular_subtype %in% c("DMG, H3 K28", "DMG, H3 wildtype") ~
        "DMG, H3 K28-altered, EZHIP",
      TRUE ~ paste0(molecular_subtype, ", EZHIP")
    )
  ) %>%
  # lesion-level refinement from the LGG lesions, for the LGG / GNT / GNG / SEGA
  # families (the LGG cohort's own call is the subtype, shown here too)
  mutate(lesion_refinement = if_else(
    lgg_specific %in% TRUE &
      (str_detect(replace_na(molecular_subtype, ""), "^(LGG|GNT|GNG|SEGA)") |
         subtype_source %in% "LGG lesion"),
    lgg_subtype, NA_character_
  )) %>%
  # what the lesion refinement could see: specific (a lesion term), none
  # detected (DNA and RNA available), or incomplete (a missing assay, so the
  # lesion is still to be determined, not absent)
  mutate(
    in_lgg_families = str_detect(replace_na(molecular_subtype, ""),
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
  # the methylation class carries its lesion term from the LGG lesion module, so
  # molecular_subtype is the lesion-level label, as in the legacy LGAT module:
  # - a generic MAPK bucket (PA_MID "LGG, MAPK", PA_INF "LGG, BRAF/MAPK", GG
  #   "GNG, MAPK") is replaced by the lesion term (LGG, KIAA1549-BRAF)
  # - PA ST gets the lesion term appended (GNG, PA ST, BRAF V600E): the class is
  #   a histology with several MAPK drivers
  # - an entity that is defined as such (DNET, SEGA, DGONC, DIG/DIA, papillary
  #   and myxoid GNT), a class that already is a lesion term (FGFR, IDH,
  #   MYB/MYBL1, PXA BRAF V600E) and the DLGNT classes are left alone; their
  #   lesions are in lesion_refinement and subtype_alterations_confirmed
  # without a lesion term the label stays at the class and
  # lesion_refinement_status says why
  mutate(
    methyl_term = str_remove(replace_na(molecular_subtype, ""), "^[^,]*,? ?"),
    methyl_prefix = str_extract(replace_na(molecular_subtype, ""), "^[^,]*"),
    refinement_term = str_remove(replace_na(lesion_refinement, ""), "^[^,]*, "),
    class_refined = subtype_source %in% "methylation" &
      !is.na(lesion_refinement) & refinement_term != "" &
      methyl_prefix %in% c("LGG", "GNG", "GNT", "SEGA") &
      !str_detect(replace_na(methyl_class, ""), "DLGNT"),
    molecular_subtype = case_when(
      class_refined & methyl_term %in% c("MAPK", "BRAF/MAPK") ~
        paste0(methyl_prefix, ", ", refinement_term),
      class_refined &
        methyl_term %in% "PA ST" &
        !str_detect(methyl_term, fixed(refinement_term)) ~
        paste0(molecular_subtype, ", ", refinement_term),
      TRUE ~ molecular_subtype
    )
  ) %>%
  select(-methyl_term, -methyl_prefix, -refinement_term, -class_refined) %>%
  # a ganglioglioma diagnosis keeps the GNG prefix for the pilocytic-type and
  # MAPK labels (KIAA1549-BRAF, BRAF V600E, other MAPK, RTK, FGFR, NF1, ...), as
  # the LGG lesion module does; IDH, MYB/MYBL1 and PXA labels are entities of their
  # own and keep LGG
  mutate(molecular_subtype = if_else(
    subtype_source %in% c("methylation", "lesion fallback") &
      pathology_diagnosis %in% "Ganglioglioma" &
      str_detect(replace_na(molecular_subtype, ""), "^LGG, ") &
      !str_detect(replace_na(molecular_subtype, ""), "^LGG, (IDH|MYB/MYBL1|PXA)"),
    str_replace(molecular_subtype, "^LGG, ", "GNG, "), molecular_subtype)) %>%
  # umbrella group (e.g. the pediatric H3-wildtype HGG classes), unless an H3
  # variant was found, which contradicts the wildtype grouping
  mutate(subtype_group = map_chr(molecular_subtype, umbrella_group)) %>%
  mutate(
    h3_group_conflict = !is.na(subtype_group) & h3_status %in% "mutant",
    subtype_group = case_when(
      h3_group_conflict ~ molecular_subtype,
      !is.na(subtype_group) ~ subtype_group,
      TRUE ~ molecular_subtype
    ),
    needs_review = is.na(manual_subtype) &
      (dx_concordance %in% "discordant" | methyl_conflict |
         h3_group_conflict | h3_variant_conflict %in% TRUE |
         mb_concordance %in% "discordant" | !is.na(lesion_contradiction) |
         !is.na(free_text_conflict) | !is.na(solid_review_flag))
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
      !str_detect(replace_na(molecular_subtype, ""), "EGFR") ~ NA_character_,
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
                  subtype_source %in% c("lesion fallback", "lesion overrides methylation",
                              "MYCN copy number", "MYCN expression",
                              "EWSR1 fusion"),
                lesion_evidence, NA_character_),
        if_else(h3_status %in% "mutant" & str_detect(replace_na(molecular_subtype, ""),
                                                     "-mutant"), h3_variant, NA_character_),
        egfr_evidence,
        if_else(str_detect(replace_na(molecular_subtype, ""), "EZHIP"),
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
           is.na(manual_subtype) & !is.na(free_text_conflict),
           # methylation and genomics agree with each other and the diagnosis
           # differs: the molecular data may not belong to the clinical record
           is.na(manual_subtype) & dx_concordance %in% "discordant" &
             lesion_status %in% "confirmed" & subtype_source %in% "methylation"),
      function(...) {
        flags <- c("diagnosis discordant", "classifier versions disagree",
                   "lesion contradicts methylation", "H3 variant outside DMG/DHG",
                   "H3 variant in wildtype group", "expression classifier discordant",
                   "free-text diagnosis differs from methylation",
                   "methylation and genomics agree, diagnosis differs")
        hit <- flags[c(...)]
        if (length(hit) == 0) NA_character_ else paste(hit, collapse = ";")
      }
    ),
    # evidence flags of the solid tumors (06-subtype-solid-tumors.R)
    subtype_review_flag = case_when(
      !is.na(solid_review_flag) & is.na(subtype_review_flag) ~ solid_review_flag,
      !is.na(solid_review_flag) ~ paste(subtype_review_flag, solid_review_flag,
                                        sep = ";"),
      TRUE ~ subtype_review_flag
    ),
    Notes = subtype_notes(subtype_source, molecular_subtype)
  ) %>%
  # tumors whose methylation conflicts with the diagnosis: `likely` where the
  # genomics supports the diagnosis (the methylation array is probably from
  # another specimen: collect and remove from the cohort), `possible` where it
  # is unresolved (input/possible_methylation_mixups.tsv)
  left_join(reviewed %>% transmute(sample_id,
                                   diagnosis_review = paste0(decision, " (", source, "): ", note)),
            by = "sample_id") %>%
  mutate(
    # a reviewed tumor has no review flag left (the review saw the flags of the
    # tumor); a methylation mismatch is listed in input/possible_methylation_mixups.tsv
    subtype_review_flag = if_else(is.na(diagnosis_review), subtype_review_flag, NA_character_),
    needs_review = needs_review & is.na(diagnosis_review)
  ) %>%
  left_join(mismatch %>% select(sample_id, methylation_mismatch = mismatch_call,
                                methylation_mismatch_evidence = genomics_support),
            by = "sample_id") %>%
  mutate(
    methylation_mismatch = if_else(methylation_mismatch %in% c("likely", "possible"),
                                   methylation_mismatch, NA_character_),
    methylation_mismatch_evidence = if_else(is.na(methylation_mismatch), NA_character_,
                                            methylation_mismatch_evidence),
    subtype_review_flag = case_when(
      is.na(methylation_mismatch) ~ subtype_review_flag,
      is.na(subtype_review_flag) ~ paste("methylation", methylation_mismatch, "mismatch"),
      TRUE ~ paste0(subtype_review_flag, ";methylation ", methylation_mismatch, " mismatch")
    )
  ) %>%
  select(match_id, Kids_First_Participant_ID, sample_id,
         molecular_subtype, subtype_group, lesion_refinement,
         lesion_refinement_status, subtype_source,
         subtype_alterations_confirmed, alteration_confirmed, subtype_review_flag, Notes,
         methylation_mismatch, methylation_mismatch_evidence, diagnosis_review,
         automated_subtype, manual_conflict, manual_source, manual_note,
         evidence_level, lesion_status, lesion_evidence, lesion_contradiction,
         methyl_class,
         mb_expression_group, mb_concordance, mb_shh_subtype,
         smarcb1_focal_status, smarcb1_expression,
         methyl_score, methyl_conflict, pathology_diagnosis, dx_group,
         dx_concordance, needs_review, free_text_conflict,
         mycn_cn_status, mycn_tpm, mycn_text_status, starts_with("tp53_"),
         ends_with("_mutation_status"), ends_with("_variants"),
         ends_with("_amplification"), ends_with("_activating_alteration"),
         ends_with("_tpm"), ends_with("_z"), ends_with("_expression_status"),
         starts_with("h3_"), starts_with("h3k27_")) %>%
  arrange(match_id) %>%
  # the DMG profile lesions (ACVR1, PPM1D, PIK3CA) of a DMG in one column
  mutate(dmg_profile_alterations = pmap_chr(
    list(acvr1_activating_variants, ppm1d_activating_variants, pik3ca_activating_variants),
    function(a, p, k) {
      x <- c(if (!is.na(a)) paste0("ACVR1 ", a), if (!is.na(p)) paste0("PPM1D ", p),
             if (!is.na(k)) paste0("PIK3CA ", k))
      if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
    }))

stopifnot(!anyDuplicated(compiled$match_id))

# the tumors with a likely methylation mismatch and their methylation arrays, to
# remove from the cohort. Built from the mismatch calls, not from the subtype table:
# a tumor with no other evidence (no DNA, no free text rule) has no subtype row but
# its arrays are still to be removed
methyl_specimens <- hist_base %>%
  filter(experimental_strategy == "Methylation") %>%
  distinct(sample_id, match_id) %>%
  distinct(sample_id, .keep_all = TRUE)
likely_mismatch <- mismatch %>%
  filter(mismatch_call == "likely") %>%
  left_join(compiled %>%
              select(sample_id, match_id, pathology_diagnosis, molecular_subtype,
                     methylation_mismatch_evidence),
            by = "sample_id") %>%
  # a tumor without a subtype row: the match_id from the methylation specimens, the
  # diagnosis and the evidence from the mismatch call
  left_join(methyl_specimens %>% rename(methyl_match_id = match_id), by = "sample_id") %>%
  mutate(match_id = coalesce(match_id, methyl_match_id),
         pathology_diagnosis = coalesce(pathology_diagnosis, diagnosis),
         methylation_mismatch_evidence = coalesce(methylation_mismatch_evidence,
                                                  conflicting_evidence)) %>%
  select(sample_id, match_id, pathology_diagnosis, molecular_subtype,
         methylation_mismatch_evidence, methyl_class, methyl_score, array_calls,
         methylation_biospecimen_ids)
stopifnot(nrow(likely_mismatch) == sum(mismatch$mismatch_call == "likely"),
          !anyDuplicated(likely_mismatch$sample_id))
write_tsv(likely_mismatch, file.path(results_dir, "methylation_likely_mismatch.tsv"))

# diffuse leptomeningeal glioneuronal tumors: a KIAA1549::BRAF fusion (or another
# MAPK alteration) with a 1p deletion is the lesion pair of the tumor
# (ploidy-aware 1p, see 02-prepare-lesion-tables.R). The subtype keeps the
# methylation or LGG lesion call; this lists the lesions.
dlgnt_arm <- read_tsv(file.path(root_dir, "scratch", "molecular-subtyping", "arm_values.tsv"),
                      show_col_types = FALSE) %>%
  select(Kids_First_Biospecimen_ID, `1p`) %>%
  inner_join(hist_base %>% distinct(Kids_First_Biospecimen_ID, match_id),
             by = "Kids_First_Biospecimen_ID") %>%
  group_by(match_id) %>%
  summarise(arm_1p = min(`1p`, na.rm = TRUE), .groups = "drop") %>%
  mutate(arm_1p = if_else(is.finite(arm_1p), arm_1p, NA_real_))
dlgnt <- hist_base %>%
  filter(sample_type == "Tumor", !is.na(match_id),
         str_detect(replace_na(pathology_free_text_diagnosis, ""),
                    regex("leptomeningeal glioneuronal", ignore_case = TRUE)) |
           (sample_type == "Tumor" & !is.na(match_id) &
              str_detect(replace_na(pathology_diagnosis, ""),
                         regex("leptomeningeal glioneuronal", ignore_case = TRUE)))) %>%
  distinct(match_id, sample_id) %>%
  left_join(compiled %>% select(match_id, molecular_subtype, lesion_evidence,
                                subtype_alterations_confirmed), by = "match_id") %>%
  left_join(dlgnt_arm, by = "match_id") %>%
  mutate(braf_fusion = str_detect(replace_na(subtype_alterations_confirmed, ""), "KIAA1549--BRAF") |
           str_detect(molecular_subtype, "KIAA1549-BRAF"),
         deletion_1p = arm_1p <= -0.5)
write_tsv(dlgnt, file.path(results_dir, "dlgnt_lesions.tsv"))

compiled <- add_cell_line_comparison(compiled, root_dir, results_dir)

write_tsv(compiled, file.path(results_dir, "molecular_subtype.tsv"))

message("Tumors with a subtype: ",
        sum(!is.na(compiled$molecular_subtype)), " of ", nrow(compiled))
print(count(compiled, subtype_source))
message("With lesion refinement: ", sum(!is.na(compiled$lesion_refinement)))
message("\nNotes")
print(count(compiled, Notes, sort = TRUE))
message("\nReview flags")
print(count(compiled, subtype_review_flag, sort = TRUE))

### Provenance of the run ---------------------------------------------------------
# the input files the results were made on, so that a rerun on other data is visible
# (utils/check-subtyping.R compares them with the files that are there)
provenance_files <- c("data/histologies-base.tsv", "data/histologies.tsv",
                      "data/fusion-putative-oncogenic.tsv",
                      "analyses/molecular-subtyping/input/bti_methylation_classifications.tsv",
                      "analyses/molecular-subtyping/input/subtype_lesion_rules.tsv")
provenance <- tibble(file = provenance_files,
                     md5 = unname(tools::md5sum(file.path(root_dir, provenance_files))),
                     release = basename(Sys.readlink(file.path(root_dir, "data", "histologies-base.tsv")) %>%
                                          dirname()))
write_tsv(provenance, file.path(results_dir, "run_provenance.tsv"))
