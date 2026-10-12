################################################################################
# 04-lgg-lesions.R
#
# LGG lesion-level subtypes (KIAA1549-BRAF, BRAF V600E, other MAPK, RTK, FGFR,
# NF1, IDH, H3, MYB/MYBL1, CDKN2A/B), ported from molecular-subtyping-LGAT so
# that this module computes the evidence itself.
#
# Components are defined in input/lgg_component_rules.tsv. As in the LGAT module:
#   - fusions come from fusion_summary_lgg_hgg_foi.tsv (in-frame, kinase domain
#     retained), and a BRAF / RAF1 fusion is "other MAPK" only when the tumor has
#     no KIAA1549--BRAF
#   - a BRAF variant is "other MAPK" only when it is not already a BRAF V600E
#     (hotspot 599/600 or kinase domain) variant
#   - NF1-germline is a cancer predisposition of NF-1
#   - FGFR also covers a FGFR1 kinase-domain tandem duplication or a copy-number
#     segment within FGFR1; CDKN2A/B is copy number <= 1
# Components are joined in the LGAT order into one label ("LGG, KIAA1549-BRAF,
# NF1-germline"), prefixed LGG / GNT / GNG / SEGA by the pathology diagnosis.
# `lgg_cohort` marks the tumors the LGAT module subtyped (its diagnosis strings);
# the lesions are computed for every PBTA tumor.
#
# Usage: Rscript --vanilla 04-lgg-lesions.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-CNS")
source(file.path(module_dir, "utils", "duckdb-connect.R"))
input_dir <- file.path(module_dir, "input")
results_dir <- file.path(module_dir, "results")
scratch_dir <- file.path(root_dir, "scratch", "molecular-subtyping-CNS")

components <- read_tsv(file.path(input_dir, "lgg_component_rules.tsv"),
                       show_col_types = FALSE, na = character())

### Tumors and assays ----------------------------------------------------------

hist <- read_tsv(file.path(root_dir, "data", "histologies-base.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  filter(cohort == "PBTA", sample_type == "Tumor", !is.na(match_id)) %>%
  select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID, sample_id,
         match_id, experimental_strategy, RNA_library, cancer_predispositions,
         pathology_diagnosis, pathology_free_text_diagnosis)

# DNA and RNA biospecimens, as in the LGAT module
dna_ids <- hist %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Sequencing"),
         is.na(RNA_library)) %>%
  distinct(match_id, id = Kids_First_Biospecimen_ID)
rna_ids <- hist %>%
  filter(experimental_strategy == "RNA-Seq" |
           (experimental_strategy == "Targeted Sequencing" & !is.na(RNA_library))) %>%
  distinct(match_id, id = Kids_First_Biospecimen_ID)

tumors <- hist %>%
  group_by(match_id) %>%
  summarise(
    Kids_First_Participant_ID = first(Kids_First_Participant_ID),
    sample_id = first(sample_id),
    pathology_diagnosis = paste(sort(unique(na.omit(pathology_diagnosis))), collapse = ";"),
    pathology_free_text_diagnosis = paste(sort(unique(na.omit(pathology_free_text_diagnosis))),
                                          collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(has_dna = match_id %in% dna_ids$match_id,
         has_rna = match_id %in% rna_ids$match_id)

### LGG cohort (the LGAT module's diagnosis strings) --------------------------
dx_strings <- jsonlite::fromJSON(file.path(input_dir, "lgg_path_dx_strings.json"))
in_cohort <- function(dx) {
  dx <- str_to_lower(replace_na(dx, ""))
  str_detect(dx, paste(dx_strings$include_path_dx, collapse = "|")) &
    !str_detect(dx, paste(dx_strings$exclude_path_dx, collapse = "|"))
}

### Evidence -------------------------------------------------------------------

con <- connect_duckdb()
snv <- dbGetQuery(con, sprintf(
  "SELECT * FROM read_parquet('%s')",
  file.path(scratch_dir, "snv_lesion_genes.parquet")
)) %>%
  as_tibble() %>%
  # synonymous calls are removed, as in the LGAT module
  filter(!is.na(HGVSp_Short), !str_detect(HGVSp_Short, "=$"))

cdkn <- dbGetQuery(con, sprintf(
  "SELECT * FROM read_parquet('%s')",
  file.path(scratch_dir, "cdkn2_copy_number.parquet")
)) %>% as_tibble()
dbDisconnect(con, shutdown = TRUE)

fgfr_dup <- read_tsv(file.path(scratch_dir, "fgfr1_dup.tsv"), show_col_types = FALSE)

# fusion summary: one column per fusion, 1 when present
fusion_long <- read_tsv(
  file.path(root_dir, "data", "fusion_summary_lgg_hgg_foi.tsv"),
  show_col_types = FALSE
) %>%
  pivot_longer(-Kids_First_Biospecimen_ID, names_to = "fusion", values_to = "n") %>%
  filter(n > 0) %>%
  inner_join(rna_ids, by = c("Kids_First_Biospecimen_ID" = "id"))

# TRUE for every row when the pattern is empty
matches <- function(x, pattern) {
  if (pattern == "") rep(TRUE, length(x)) else str_detect(x, pattern)
}

# one (match_id, component, evidence) row per piece of evidence
snv_hits <- components %>%
  filter(lesion_type == "snv") %>%
  mutate(rule = row_number()) %>%
  pmap_dfr(function(component, gene, pattern, variant_class, domain, rule, ...) {
    snv %>%
      filter(str_detect(Hugo_Symbol, paste0("^(", gene, ")$")),
             str_detect(HGVSp_Short, pattern),
             matches(Variant_Classification, variant_class),
             matches(replace_na(DOMAINS, ""), domain)) %>%
      transmute(id = Kids_First_Biospecimen_ID, component,
                evidence = paste(Hugo_Symbol, HGVSp_Short),
                is_v600e_rule = component == "BRAF V600E")
  }) %>%
  inner_join(dna_ids, by = "id", relationship = "many-to-many")

# BRAF variants that are not BRAF V600E become "other MAPK"
braf_v600e_ids <- snv_hits %>% filter(component == "BRAF V600E") %>%
  distinct(id, evidence)
braf_other <- snv %>%
  filter(Hugo_Symbol == "BRAF") %>%
  anti_join(braf_v600e_ids %>% mutate(evidence = evidence),
            by = c("Kids_First_Biospecimen_ID" = "id")) %>%
  transmute(id = Kids_First_Biospecimen_ID, component = "other MAPK",
            evidence = paste(Hugo_Symbol, HGVSp_Short)) %>%
  inner_join(dna_ids, by = "id", relationship = "many-to-many")

fusion_hits <- components %>%
  filter(lesion_type == "fusion") %>%
  pmap_dfr(function(component, pattern, ...) {
    fusion_long %>%
      filter(str_detect(fusion, pattern)) %>%
      transmute(match_id, component, evidence = paste0(fusion, " fusion"))
  })

# other MAPK fusions do not count when the tumor has KIAA1549--BRAF
canonical <- fusion_hits %>% filter(component == "KIAA1549-BRAF") %>%
  distinct(match_id)
fusion_hits <- fusion_hits %>%
  filter(!(component == "other MAPK" & match_id %in% canonical$match_id))

nf1_germline <- hist %>%
  filter(experimental_strategy != "RNA-Seq",
         str_detect(replace_na(cancer_predispositions, ""), "NF-1")) %>%
  distinct(match_id) %>%
  transmute(match_id, component = "NF1-germline",
            evidence = "NF-1 cancer predisposition")

cnv_hits <- bind_rows(
  cdkn %>% filter(copy_number <= 1) %>%
    transmute(id = Kids_First_Biospecimen_ID, component = "CDKN2A/B",
              evidence = paste0(gene_symbol, " copy number ", copy_number)) %>%
    inner_join(dna_ids, by = "id", relationship = "many-to-many"),
  fgfr_dup %>%
    transmute(id = Kids_First_Biospecimen_ID, component = "FGFR",
              evidence = paste0("FGFR1 ", fgfr1_dup_source)) %>%
    inner_join(dna_ids, by = "id", relationship = "many-to-many")
)

hits <- bind_rows(
  snv_hits %>% select(match_id, component, evidence),
  braf_other %>% select(match_id, component, evidence),
  fusion_hits,
  nf1_germline,
  cnv_hits %>% select(match_id, component, evidence)
) %>%
  # NF1 rows from the SNV rule are somatic
  mutate(component = recode(component, "NF1-somatic" = "NF1-somatic")) %>%
  distinct()

### Labels (components in the LGAT order) --------------------------------------

# the LGAT module lists RNA-derived components first and DNA-derived ones next,
# each in a fixed order, and keeps the first occurrence of a component
rna_order <- c("KIAA1549-BRAF", "other MAPK", "RTK", "FGFR", "MYB/MYBL1")
dna_order <- c("NF1-somatic", "NF1-germline", "BRAF V600E", "other MAPK", "RTK",
               "FGFR", "IDH", "H3", "CDKN2A/B")
rank_levels <- c(paste0("rna|", rna_order), paste0("dna|", dna_order))

labels <- hits %>%
  mutate(assay = if_else(component %in% rna_order & str_detect(evidence, "fusion$"),
                         "rna", "dna"),
         rank = match(paste0(assay, "|", component), rank_levels)) %>%
  filter(!is.na(rank)) %>%
  group_by(match_id) %>%
  arrange(rank, .by_group = TRUE) %>%
  summarise(
    components = paste(unique(component), collapse = ", "),
    lgg_alterations = paste(sort(unique(evidence)), collapse = "; "),
    .groups = "drop"
  )

out <- tumors %>%
  left_join(labels, by = "match_id") %>%
  mutate(
    lgg_lesion_status = case_when(
      !is.na(components) ~ "specific",
      has_dna & has_rna ~ "none detected",
      TRUE ~ "incomplete"
    ),
    prefix = case_when(
      pathology_diagnosis == "Glial-neuronal tumor NOS" ~ "GNT",
      str_detect(str_to_lower(replace_na(pathology_free_text_diagnosis, "")),
                 "desmoplastic infantile|glioneuronal") &
        pathology_diagnosis == "Low-grade glioma/astrocytoma (WHO grade I/II)" ~ "GNT",
      pathology_diagnosis == "Ganglioglioma" ~ "GNG",
      pathology_diagnosis == "Subependymal Giant Cell Astrocytoma (SEGA)" ~ "SEGA",
      TRUE ~ "LGG"
    ),
    lgg_lesion_label = if_else(lgg_lesion_status == "specific",
                               paste0(prefix, ", ", components), NA_character_),
    lgg_cohort = in_cohort(pathology_diagnosis)
  ) %>%
  select(match_id, Kids_First_Participant_ID, sample_id, has_dna, has_rna,
         lgg_cohort, lgg_lesion_status, lgg_lesion_label, lgg_alterations)

write_tsv(out, file.path(results_dir, "lgg_lesions.tsv"))

message("Tumors: ", nrow(out))
print(count(out, lgg_lesion_status))
print(out %>% filter(lgg_cohort, lgg_lesion_status == "specific") %>%
        count(lgg_lesion_label, sort = TRUE) %>% head(15))
