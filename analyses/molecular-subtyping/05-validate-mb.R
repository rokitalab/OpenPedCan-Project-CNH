################################################################################
# 05-validate-mb.R
#
# Medulloblastoma. Three steps:
#   1. the RNA-seq expression classifier (medulloPackage) on the tumors with a
#      medulloblastoma diagnosis;
#   2. the SHH subtype (alpha, beta, gamma, delta): the methylation subclass
#      (DKFZ v12.8 MB_SHH_1 to 4 at >= 0.7, then NIH v3.1 MB SHH 1 to 4 at
#      >= 0.9) where there is one, otherwise age and molecular criteria (copy
#      number, SNV, germline, expression);
#   3. comparison of the methylation group with the expression classifier group.
#
# Methylation sets the group (03-validate-subtypes.R, with lesion confirmation
# for WNT, SHH and Group 3). The expression classifier is an independent check,
# and the subtype for tumors with no methylation call. Disagreements are flagged.
#
# Usage: Rscript --vanilla 05-validate-mb.R
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DBI)
  # attached, not only called with ::, because classify() loads its gene sets
  # with utils::data("medulloSetsUp"), which needs the package on the search path
  library(medulloPackage)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping")
results_dir <- file.path(module_dir, "results")
input_dir <- file.path(module_dir, "input")
data_dir <- file.path(root_dir, "data")
source(file.path(module_dir, "utils", "duckdb-connect.R"))
source(file.path(module_dir, "utils", "mismatch-arrays.R"))
source(file.path(module_dir, "utils", "release-tables.R"))

hist <- read_tsv(file.path(data_dir, "histologies-base.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  select(-any_of("molecular_subtype")) %>%
  mutate(age_at_diagnosis_years = age_at_diagnosis_days / 365.25)

calls <- read_tsv(file.path(results_dir, "methylation_calls.tsv"),
                  show_col_types = FALSE) %>%
  exclude_mismatched_arrays(results_dir)

### Cohort ---------------------------------------------------------------------
# Medulloblastoma diagnosis, or "Other" with medulloblastoma in the free text
mb_dx <- jsonlite::fromJSON(file.path(input_dir, "path_dx_strings.json"))$mb
mb_biospecimens <- hist %>%
  filter(pathology_diagnosis %in% mb_dx$exact_path_dx |
           (pathology_diagnosis == "Other" &
              str_detect(str_to_lower(pathology_free_text_diagnosis),
                         paste0(mb_dx$include_free_text, collapse = "|")))) %>%
  filter(!is.na(match_id))

### Expression classifier ------------------------------------------------------
# log2(TPM + 1) of the RNA-seq biospecimens of the cohort
expr <- readRDS(file.path(data_dir, "gene-expression-rsem-tpm-collapsed.rds"))
mb_rna_ids <- mb_biospecimens %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  pull(Kids_First_Biospecimen_ID) %>%
  intersect(colnames(expr))
mb_classified <- medulloPackage::classify(log2(expr[, mb_rna_ids] + 1)) %>%
  arrange(sample)
write_tsv(mb_classified, file.path(results_dir, "mb_expression_classifier.tsv"))

# one call per RNA biospecimen, collapsed to the tumor
expression <- mb_classified %>%
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

### Methylation group and lesion confirmation from the validator ---------------
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


### SHH subtype ----------------------------------------------------------------
# SHH tumors: methylation group MB, SHH, or no methylation group and the
# expression classifier says SHH
shh_match_ids <- methylation %>%
  full_join(expression, by = "match_id") %>%
  mutate(group = coalesce(mb_methylation_group, mb_expression_group)) %>%
  filter(group == "MB, SHH") %>%
  pull(match_id)

mb_shh <- mb_biospecimens %>%
  filter(match_id %in% shh_match_ids) %>%
  distinct(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, match_id,
           sample_id) %>%
  left_join(select(methylation, match_id,
                   molecular_subtype_methyl = mb_methylation_group),
            by = "match_id") %>%
  mutate(molecular_subtype = "MB, SHH")

# SHH subclass (SHH_1 beta, SHH_2 gamma, SHH_3 alpha, SHH_4 delta) from the
# DKFZ v12.8 subclass when the mean score of the tumor's arrays is >= 0.7 (lower
# than the 0.8 used for the methylation group, as in the SHH subtyping of the
# legacy module), otherwise from the NIH v3.1 subclass (MB SHH 1 to 4, same
# numbering) when its mean score is >= 0.9. A tumor with conflicting subclasses
# in a classifier has none from it.
shh_number <- c(`1` = "SHH_beta", `2` = "SHH_gamma", `3` = "SHH_alpha",
                `4` = "SHH_delta")

bti_shh <- read_tsv(
  file.path(input_dir, "bti_methylation_classifications.tsv"),
  col_select = c(Kids_First_Biospecimen_ID = Bioassay_ID,
                 dkfz_subclass = `dkfz_v12.8_subclass`,
                 dkfz_score = `dkfz_v12.8_subclass_maxscore`,
                 nih_subclass = `NIH_v3.1_subclass`,
                 nih_score = `NIH_v3.1_subclass_meanscore`),
  show_col_types = FALSE) %>%
  inner_join(select(hist, Kids_First_Biospecimen_ID, match_id),
             by = "Kids_First_Biospecimen_ID") %>%
  filter(match_id %in% shh_match_ids)

shh_subclass <- function(df, subclass, score, pattern, cutoff, source) {
  df %>%
    transmute(match_id, Kids_First_Biospecimen_ID,
              subclass = {{ subclass }}, score = {{ score }}) %>%
    filter(!is.na(subclass)) %>%
    group_by(match_id) %>%
    summarise(subclass = if (n_distinct(subclass) == 1) first(subclass)
                         else NA_character_,
              score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(subclass), score >= cutoff,
           str_detect(subclass, pattern)) %>%
    transmute(match_id,
              SHH_subtype = unname(shh_number[str_extract(subclass, "[1-4]$")]),
              shh_methylation_source = source)
}

shh_methylation <- bind_rows(
  shh_subclass(bti_shh, dkfz_subclass, dkfz_score, "^MB_SHH_[1-4]$", 0.7,
               "dkfz_v12.8"),
  shh_subclass(bti_shh, nih_subclass, nih_score, "^MB SHH [1-4]$", 0.9,
               "nih_v3.1")
) %>%
  group_by(match_id) %>%
  slice(1) %>%
  ungroup()

# set file paths
cn_arm_file <- file.path(data_dir, "cnv-consensus-gistic.zip")
germline_file <- file.path(input_dir, "pbta-merged-plp-variants-cbtn-plus-mb-x01.tsv")
germline_sv_file <- file.path(input_dir, "pbta_germline_svs.tsv")
tp53_alterations_file <- file.path(root_dir, "analyses", "tp53_nf1_score",
                                   "results", "tp53_altered_status.tsv")

# The gene copy number and SNV tables are large (tens of millions of rows), so
# only the rows of the SHH biospecimens and genes of interest are queried, from
# the parquet copies in data/parquet (scripts/download-parquet.sh) or else from
# the release files.
con <- connect_duckdb()
shh_ids_sql <- paste(sprintf("'%s'", mb_shh$Kids_First_Biospecimen_ID),
                     collapse = ", ")

# Load gene CN data and subset for GOIs
cn_gene <- dbGetQuery(con, sprintf(
  "SELECT biospecimen_id, gene_symbol, lower(status) AS status FROM %s
   WHERE biospecimen_id IN (%s)
     AND gene_symbol IN ('MYCN', 'GLI2', 'CCND2', 'PTEN', 'PPM1D', 'MDM4', 'TP53')
     AND lower(status) IN ('amplification', 'gain', 'loss', 'deep deletion')",
  release_source("consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only",
                 "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz",
                 root_dir),
  shh_ids_sql)) %>%
  as_tibble()

# load the SNV tables (consensus and tumor-only) and subset for GOIs
keep_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Start_Position", "End_Position",
               "HGVSp_Short", "HGVSc", "Consequence", "Variant_Classification")

snv_query <- function(name, tsv_file) {
  dbGetQuery(con, sprintf(
    "SELECT %s FROM %s
     WHERE Tumor_Sample_Barcode IN (%s)
       AND (Hugo_Symbol IN ('CREBBP', 'DDX3X', 'KMT2D', 'PTCH1', 'SMO', 'SUFU', 'TERT')
            -- for U1 mutations, keep only those at 3rd and 5th nucleotide
            OR (regexp_matches(Hugo_Symbol, 'RN*U1-*')
                AND HGVSc IN ('n.3A>G', 'n.5A>G')))",
    paste(keep_cols, collapse = ", "),
    release_source(name, tsv_file, root_dir), shh_ids_sql)) %>%
    as_tibble() %>%
    mutate(across(c(Start_Position, End_Position), as.character))
}
mutations <- bind_rows(
  snv_query("snv-consensus-plus-hotspots", "snv-consensus-plus-hotspots.maf.tsv.gz"),
  snv_query("snv-mutect2-tumor-only-plus-hotspots",
            "snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz"))
dbDisconnect(con, shutdown = TRUE)

# get list of U1 genes in this dataset
u1_genes <- mutations %>%
  filter(grepl("RN*U1-*", Hugo_Symbol)) %>%
  pull(Hugo_Symbol) %>%
  unique()

# Load CN chromosome arm file
broad_CNA <- readr::read_tsv(unz(cn_arm_file, "cnv-consensus-gistic/broad_values_by_arm.txt")) %>%
  column_to_rownames("Chromosome Arm")

# Define function for filling in alterations data frame
fill_df <- function(sample, ref_df, col_name, included_samples = NULL, default = 0){
  if(is.na(sample)) { # no results for this sample
    return(NA)
  } else if(!is.null(included_samples) && !sample %in% included_samples) { # sample is not present in included samples, return NA
    return(NA)
  } else if(!sample %in% rownames(ref_df)){ # sample is not present in the reference data set, return default
    return(default)
  } else {
    value <- ref_df[sample, col_name]
    return(value)
  }
}

broad_CNA_fill_df <- function(x, broad_CNA_data, arm, loss_gain = c("loss", "gain")){
  biospecimen_DNA <- x[['Kids_First_Biospecimen_ID']]
  # no DNA results for this sample
  if(is.na(biospecimen_DNA)){
    return(NA)
  }
  CNA_value <- broad_CNA_data[arm, biospecimen_DNA]
  if(is.na(CNA_value) || is.null(CNA_value)){
    return(NA)
  } else if(CNA_value < 0 && loss_gain == "loss") {
    return("1")
  } else if(CNA_value > 0 && loss_gain == "gain"){
    return(1)
  } else {
    return(0)
  }
}

mutation_fill_df <- function(sample, ref_df, col_name, default = 0){
  ref_df <-  ref_df %>% as_tibble()
  if(is.na(sample)) { # no results for this sample
    return(NA)
  } else if(!sample %in% ref_df$Tumor_Sample_Barcode){ # sample is not present in the reference data set, return default
    return(default)
  } else {
    value <- ref_df[ref_df$Tumor_Sample_Barcode == sample, col_name]
    return(value)
  }
}

# build molecular alteration data frame 
mb_shh_data <- mb_shh %>%
  dplyr::select(Kids_First_Biospecimen_ID, Kids_First_Participant_ID, match_id)

# add arm gain/loss info
for(arm_of_interest in c("2q", "2p", "9q", "9p", "10q", "14q", "17p")) {
  mb_shh_data[paste0(arm_of_interest, "_loss")] <- apply(mb_shh_data, MARGIN = 1, FUN = function(x) 
    broad_CNA_fill_df(x = x, broad_CNA_data = broad_CNA, arm = arm_of_interest, loss_gain = "loss"))
  mb_shh_data[paste0(arm_of_interest, "_gain")] <- apply(mb_shh_data, MARGIN = 1, FUN = function(x) 
    broad_CNA_fill_df(x = x, broad_CNA_data = broad_CNA, arm = arm_of_interest, loss_gain = "gain"))
}

# Add CN data
for (goi in c("MYCN", "GLI2", "CCND2", "PTEN", "PPM1D", "MDM4", "TP53")) {
  cn_gene_filtered <- cn_gene %>%
    filter(gene_symbol == goi) %>%
    column_to_rownames("biospecimen_id")
  mb_shh_data[,paste0("consensus_CN_", goi)] <- sapply(mb_shh_data$Kids_First_Biospecimen_ID, 
                                                      function(x) fill_df(sample = x, ref_df = cn_gene_filtered, 
                                                                          col_name = "status", default = NA_character_))
}

# Loop through GOIs to add SNV/INDEL consequences to data frame
for (gene in c("CREBBP", "DDX3X", "KMT2D", "PTCH1", "SMO", "SUFU", "TERT", "TP53", u1_genes)){
  
  if (gene == "TERT"){
    
    mutation_gene <- mutations %>% filter(Hugo_Symbol == "TERT",
                                          grepl("5'Flank", Variant_Classification),
                                          Start_Position %in% c("1295113","1295135"),
                                          End_Position %in% c("1295113","1295135")) %>%
      group_by(Tumor_Sample_Barcode) %>%
      summarise(Consequences = str_c(Consequence, collapse = ";"))
    
  } 
  
  if (grepl("RN*U1-*", gene)){
    
    mutation_gene <- mutations %>% filter(Hugo_Symbol == gene & 
                                            HGVSc %in% c("n.3A>G", "n.5A>G")) %>%
      group_by(Tumor_Sample_Barcode) %>%
      summarise(Consequences = str_c(Consequence, collapse = ";"))
    
  } 
  
  else{
    
    mutation_gene <- mutations %>% filter(Hugo_Symbol == gene & grepl("frameshift|missense|stop|splice_acceptor|splice_donor", Consequence)) %>%
      group_by(Tumor_Sample_Barcode) %>%
      summarise(Consequences = str_c(Consequence, collapse = ";"))
    
  }

  mb_shh_data[,paste0(gene, "_csq")] <- unlist(sapply(mb_shh_data$Kids_First_Biospecimen_ID,
                                              function(x) mutation_fill_df(sample = x, ref_df = mutation_gene, 
                                                                           col_name = 'Consequences', default = NA)))
  
}

# Collapse somatic alteration data per match ID
mb_shh_data_collapsed <- mb_shh_data %>%
  group_by(Kids_First_Participant_ID, match_id) %>%
  summarize(
    across(`2q_loss`:last_col(), ~paste0(unique(na.omit(.)), collapse = ':'), .names = "{.col}"),
    .groups = "drop") %>%
  mutate(across(everything(), ~str_remove_all(., "NA:|:NA"))) %>%
  mutate(across(everything(), ~if_else(str_detect(., "^$"), NA_character_, .))) %>%
  # add column for all U1 genes' annotation for easier selection later
  mutate(U1_mutation = as.integer(rowSums(across(matches("RN.*U1-.*"), ~ !is.na(.) & . != ""), na.rm = TRUE) > 0))

# Load germline SNV/INDEL and SV data, filter for GOIs, and merge
germline <- read_tsv(germline_file) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_normal) %>%
  left_join(hist %>% dplyr::select("Kids_First_Biospecimen_ID",
                                   "Kids_First_Participant_ID")) %>%
  dplyr::filter(Kids_First_Participant_ID %in% mb_shh_data$Kids_First_Participant_ID,
                gene_symbol_vep %in% c("ELP1", "PTCH1",
                                       "SUFU", "TP53")) %>%
  dplyr::select(Kids_First_Participant_ID, gene_symbol_vep, variant_classification_vep) %>%
  dplyr::rename(consequence = variant_classification_vep)

germline_sv <- read_tsv(germline_sv_file) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_normal) %>%
  left_join(hist %>% dplyr::select("Kids_First_Biospecimen_ID",
                                   "Kids_First_Participant_ID")) %>%
  dplyr::filter(Kids_First_Participant_ID %in% mb_shh_data$Kids_First_Participant_ID,
                Hugo_Symbol_cpg %in% c("ELP1", "PTCH1",
                                       "SUFU", "TP53")) %>%
  dplyr::select(Kids_First_Participant_ID, Hugo_Symbol_cpg, Type) %>%
  dplyr::rename(gene_symbol_vep = Hugo_Symbol_cpg,
                consequence = Type)

germline <- germline %>%
  bind_rows(germline_sv) %>%
  pivot_wider(names_from = gene_symbol_vep, values_from = consequence) %>%
  rename_with(~paste0(., "_germline_plp"), -1)

genes_of_interest = c("MYCN", "GLI2", "CCND2", "PTEN", "TP53")
subset_expr <- expr[genes_of_interest, mb_shh_data$Kids_First_Biospecimen_ID[mb_shh_data$Kids_First_Biospecimen_ID %in% colnames(expr)]]

# calculate max TPM per match ID and calculate z-scores
expr_df <- as.data.frame(t(subset_expr)) %>%
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::rename(MYCN_TPM = MYCN,
                GLI2_TPM = GLI2,
                CCND2_TPM = CCND2,
                PTEN_TPM = PTEN,
                TP53_TPM = TP53) %>%
  left_join(mb_shh_data) %>%
  group_by(Kids_First_Participant_ID, match_id) %>%
  summarize(MYCN_TPM = max(MYCN_TPM),
            GLI2_TPM = max(GLI2_TPM),
            CCND2_TPM = max(CCND2_TPM),
            PTEN_TPM = min(PTEN_TPM),
            TP53_TPM = min(TP53_TPM),
            .groups = "drop") %>%
  mutate(MYCN_TPM_zscore = as.vector(scale(MYCN_TPM)),
         GLI2_TPM_zscore = as.vector(scale(GLI2_TPM)),
         CCND2_TPM_zscore = as.vector(scale(CCND2_TPM)),
         PTEN_TPM_zscore = as.vector(scale(PTEN_TPM)),
         TP53_TPM_zscore = as.vector(scale(TP53_TPM))) %>%
  unique() %>%
  mutate(MYCN_exp_status = ifelse(MYCN_TPM_zscore >= 2, "MYCN high", NA_character_),
         GLI2_exp_status = ifelse(GLI2_TPM_zscore >= 2, "GLI2 high", NA_character_),
         CCND2_exp_status = ifelse(CCND2_TPM_zscore >= 2, "CCND2 high", NA_character_),
         PTEN_exp_status = ifelse(PTEN_TPM_zscore < -2, "PTEN low", NA_character_),
         TP53_exp_status = ifelse(TP53_TPM_zscore < -2, "TP53 low", NA_character_))

tp53_df <- read_tsv(tp53_alterations_file) %>%
  dplyr::rename(tp53_SV_type = SV_type,
                tp53_SV_counts = SV_counts,
                tp53_Fusions = Fusions,
                tp53_Fusion_counts = Fusion_counts,
                tp53_cn_less_ploidy = cn_less_ploidy,
                tp53_SNV_indel_counts = SNV_indel_counts,
                tp53_CNV_loss_counts = CNV_loss_counts,
                tp53_HGVSp_Short = HGVSp_Short,
                tp53_CNV_loss_evidence = CNV_loss_evidence,
                tp53_hotspot = hotspot,
                tp53_activating = activating,
                tp53_overlap_domain = overlap_domain) %>%
  dplyr::select(match_id, cancer_predispositions, starts_with("tp53"))
  

# Create the subtypes df that:
# 1) takes the SHH subclass from methylation
# 2) appends the molecular alteration data
# 3) classifies tumors without a methylation subclass into one of the four SHH
#    subtypes using molecular and age criteria
hist_mb_shh <- hist %>%
  filter(Kids_First_Biospecimen_ID %in% mb_shh$Kids_First_Biospecimen_ID) %>%
  distinct(Kids_First_Participant_ID, sample_id, match_id, age_at_diagnosis_years) %>%
  left_join(shh_methylation, by = "match_id")

mb_subtypes <- mb_shh %>%
  select(match_id, molecular_subtype, molecular_subtype_methyl) %>%
  unique()

final_subtypes <- mb_subtypes %>%
  left_join(hist_mb_shh) %>%
  left_join(mb_shh_data_collapsed) %>%
  left_join(germline) %>%
  left_join(expr_df) %>%
  left_join(tp53_df) %>%
  unique() %>%
  mutate(final_shh_subtype = case_when(
    # For SHH alpha: >3 years AND (MYCN OR GLI2 amplification OR ELP1 germline P/LP OR TP53 germline P/LP OR TPM>=2 OR chr17p loss OR 9p gain)
    SHH_subtype == "SHH_alpha" | (age_at_diagnosis_years >= 2 &
                                    is.na(SHH_subtype) &
                                    (grepl("amplification", consensus_CN_MYCN) |
                                       grepl("amplification", consensus_CN_GLI2) |
                                       grepl("amplification", consensus_CN_CCND2) |
                                       !is.na(ELP1_germline_plp) |
                                       !is.na(TP53_germline_plp) |
                                       MYCN_TPM_zscore >= 2 |
                                       GLI2_TPM_zscore >= 2 |
                                       CCND2_TPM_zscore >= 2 |
                                       TP53_TPM_zscore < -2 |
                                       tp53_hotspot > 0 | 
                                       `17p_loss` == 1 |
                                       grepl("1", `9p_gain`))) ~ "SHH alpha",
    
    # For SHH delta: (>=10 years AND (TERT OR DDX3X OR SMO mutations OR 14q-))
    SHH_subtype == "SHH_delta" |
      (age_at_diagnosis_years >= 10 & (U1_mutation == 1 | !is.na(TERT_csq) | !is.na(DDX3X_csq) | !is.na(SMO_csq) | `14q_loss` == 1)) ~ "SHH delta",
    
    # For SHH beta: <5 years AND (KMT2D mutations OR PTEN loss/deletion OR PTEN TPM<-2 OR chr2 gain)
    SHH_subtype == "SHH_beta" | (age_at_diagnosis_years < 5 &
                                   is.na(SHH_subtype) &
                                   (!is.na(KMT2D_csq) |
                                      grepl("loss|deep deletion", consensus_CN_PTEN)|
                                      PTEN_TPM_zscore < -2) |
                                   (`2p_gain` == 1 & `2q_gain` == 1)) ~ "SHH beta",
    TRUE ~ NA_character_)) %>%
  
    # Add gamma subtype for high confidence methylation samples OR
    # If not already subtyped, and age <5 and no 2p gain, then SHH gamma
  mutate(final_shh_subtype  = case_when(SHH_subtype == "SHH_gamma" | (is.na(final_shh_subtype) & 
                                                                        age_at_diagnosis_years < 5 &
                                                                        `2p_gain` == 0) ~ "SHH gamma",
                                        TRUE ~ final_shh_subtype)) %>%
  dplyr::mutate(classification_source = case_when(
    !is.na(SHH_subtype) ~ "Methylation",
    !is.na(final_shh_subtype) ~ "Genomic/Expression",
    TRUE ~ NA_character_
  )) %>%
  dplyr::mutate(tp53_status = case_when(
    # check for activating mutations
    tp53_activating == 1 ~ "TP53 activation",
    # Loss
    # germline P/LP
    !is.na(TP53_germline_plp) | 
    # other hotspot mutation  
    tp53_hotspot >= 1 | 
    # SNV + (CNV|SV)
    tp53_SNV_indel_counts >= 1 & (tp53_CNV_loss_counts >=1 | tp53_SV_counts >=1) |
    # >1 SNV/indel
    tp53_SNV_indel_counts > 1 |
    # frameshift or nonsense variant overlaps functional domain
    tp53_overlap_domain >=1 |
    # LFS + SNV/indel/CNV/SV
    grepl("Li-Fraumeni syndrome",cancer_predispositions) & (tp53_SNV_indel_counts >= 1 | tp53_CNV_loss_counts >= 1 | tp53_SV_counts >=1) ~ "TP53 loss",
    TRUE ~ "TP53 wildtype"
  )) %>%
  right_join(mb_shh) %>%
  # re-annotate molecular subtype
  mutate(molecular_subtype = case_when(!is.na(final_shh_subtype) ~ paste0("MB, ", final_shh_subtype), 
                                       TRUE ~ molecular_subtype)) %>%
  select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, match_id, molecular_subtype, molecular_subtype_methyl, final_shh_subtype, everything()) %>%
  arrange(Kids_First_Biospecimen_ID) %>%
  write_tsv(file.path(results_dir, "mb_shh_subtypes_w_molecular_data.tsv"))

final_subtypes %>%
  select("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID", "match_id", "sample_id", 
         "molecular_subtype", "molecular_subtype_methyl", "tp53_status",
         classification_source) %>%
  write_tsv(file.path(results_dir, "mb_shh_molecular_subtypes.tsv"))


print(as.data.frame(table(final_subtypes$molecular_subtype)))
print(table(final_subtypes$molecular_subtype, final_subtypes$classification_source))

# SHH subtype per tumor, for the validation below
shh <- final_subtypes %>%
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  summarise(mb_shh_subtype = first(molecular_subtype),
            mb_shh_source = first(classification_source),
            .groups = "drop")

### Validation -----------------------------------------------------------------
# a discordance with the expression classifier is flagged only when the
# classifier is confident (p < this value)
strong_expression_p <- 0.001

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
      # the expression classifier is not confident (p >= 0.001): methylation
      # decides and the tumor is not flagged (the UMAPs of the Group 3/4 and
      # SHH to Group 3 changes agree with the methylation calls)
      mb_expression_p >= strong_expression_p ~ "discordant, weak expression call",
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
message("\nDKFZ SHH class vs the SHH subtype from mb/")
print(mb %>% filter(!is.na(mb_shh_subtype),
                    str_detect(replace_na(methyl_class, ""), "SHH")) %>%
        mutate(class = toupper(gsub("[ _]+", " ", methyl_class))) %>%
        count(class, mb_shh_subtype) %>%
        pivot_wider(names_from = mb_shh_subtype, values_from = n, values_fill = 0),
      width = Inf)
