################################################################################
# 06-validate-mb.R
#
# Medulloblastoma: compare the methylation group with the RNA-seq expression
# classifier (medulloPackage, from molecular-subtyping-MB) and add the SHH
# subtype (alpha, beta, gamma, delta) from molecular-subtyping-MB.
#
# Methylation sets the group (03-validate-subtypes.R, with lesion confirmation
# for WNT, SHH and Group 3). The expression classifier is an independent check,
# and the subtype for tumors with no methylation call. Disagreements are flagged.
#
# Usage: Rscript --vanilla 06-validate-mb.R
################################################################################

suppressPackageStartupMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS")
results_dir <- file.path(module_dir, "results")
mb_dir <- file.path(root_dir, "analyses", "molecular-subtyping-MB", "results")

hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  select(Kids_First_Biospecimen_ID, match_id, sample_id, sample_type,
         pathology_diagnosis)

# expression classifier: one call per RNA biospecimen, collapsed to the tumor
expression <- readRDS(file.path(mb_dir, "mb-classified.rds"))$medulloPackage %>%
  transmute(Kids_First_Biospecimen_ID = sample,
            mb_expression_group = paste0("MB, ", best.fit),
            mb_expression_p = p.value) %>%
  inner_join(select(hist, Kids_First_Biospecimen_ID, match_id),
             by = "Kids_First_Biospecimen_ID") %>%
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(
    mb_expression_group = paste(sort(unique(mb_expression_group)), collapse = " | "),
    mb_expression_p = min(mb_expression_p),
    .groups = "drop"
  )

# SHH subtype from the MB module
shh <- read_tsv(file.path(mb_dir, "mb_shh_molecular_subtypes.tsv"),
                show_col_types = FALSE) %>%
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(mb_shh_subtype = first(molecular_subtype),
            mb_shh_source = first(classification_source),
            .groups = "drop")

# methylation group and lesion confirmation from the validator
methylation <- read_tsv(
  file.path(results_dir, "methylation_subtype_validation.tsv"),
  guess_max = 10000,
  col_types = cols(.default = col_guess(), tp53_snv_evidence = col_character(),
                   tp53_cnv_evidence = col_character(),
                   egfr_variants = col_character())
) %>%
  filter(str_detect(OPC_molecular_subtype, "^MB, ")) %>%
  # tumors without a match_id cannot be linked to the RNA classifier
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  arrange(desc(methyl_score), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(match_id, sample_id, mb_methylation_group = OPC_molecular_subtype,
         methyl_class, methyl_score, lesion_status, lesion_evidence,
         pathology_diagnosis)

mb <- full_join(methylation, expression, by = "match_id") %>%
  # keep tumors that are MB by methylation or by the expression classifier
  left_join(shh, by = "match_id") %>%
  left_join(hist %>% filter(sample_type == "Tumor", !is.na(match_id)) %>%
              distinct(match_id, sample_id_hist = sample_id), by = "match_id") %>%
  mutate(
    sample_id = coalesce(sample_id, sample_id_hist),
    # "MB, MYO" style rare variants are not an expression-classifier group
    mb_concordance = case_when(
      !is.na(mb_methylation_group) & is.na(mb_expression_group) ~
        "methylation only",
      is.na(mb_methylation_group) ~ "expression only",
      mb_methylation_group == mb_expression_group ~ "concordant",
      TRUE ~ "discordant"
    ),
    mb_final_group = coalesce(mb_methylation_group, mb_expression_group),
    mb_group_source = if_else(!is.na(mb_methylation_group), "methylation",
                              "expression classifier"),
    needs_review = mb_concordance == "discordant"
  ) %>%
  select(-sample_id_hist) %>%
  arrange(mb_concordance, match_id)

write_tsv(mb, file.path(results_dir, "mb_validation.tsv"))
mb %>%
  filter(needs_review) %>%
  write_tsv(file.path(results_dir, "mb_discordant.tsv"))

message("MB tumors: ", nrow(mb))
print(count(mb, mb_concordance))
message("\nMethylation vs expression classifier")
print(mb %>% filter(mb_concordance %in% c("concordant", "discordant")) %>%
        count(mb_methylation_group, mb_expression_group) %>%
        pivot_wider(names_from = mb_expression_group, values_from = n,
                    values_fill = 0), width = Inf)
message("\nLesion confirmation by group (methylation-called)")
print(mb %>% filter(!is.na(mb_methylation_group)) %>%
        count(mb_methylation_group, lesion_status) %>%
        pivot_wider(names_from = lesion_status, values_from = n, values_fill = 0))
message("\nDiscordant tumors")
print(mb %>% filter(needs_review) %>%
        select(sample_id, mb_methylation_group, methyl_class, methyl_score,
               mb_expression_group, mb_expression_p, lesion_status,
               lesion_evidence, pathology_diagnosis), width = Inf)
message("\nDKFZ SHH class vs the SHH subtype from molecular-subtyping-MB")
print(mb %>% filter(!is.na(mb_shh_subtype),
                    str_detect(replace_na(methyl_class, ""), "SHH")) %>%
        mutate(class = toupper(gsub("[ _]+", " ", methyl_class))) %>%
        count(class, mb_shh_subtype) %>%
        pivot_wider(names_from = mb_shh_subtype, values_from = n, values_fill = 0),
      width = Inf)
