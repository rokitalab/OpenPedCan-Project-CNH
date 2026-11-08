# Author: Ryan Corbett & Jo Lynne Rokita
# Function: Script to subtype MB SHH tumors 

# Load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set file paths
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "molecular-subtyping-MB")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")

# set file paths
hist_file <- file.path(data_dir, "histologies-base.tsv")
mb_file <- file.path(results_dir, "MB_molecular_subtype.tsv")
cn_file <- file.path(data_dir, "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only_autosomes.tsv.gz")
maf_file <- file.path(data_dir, "snv-consensus-plus-hotspots.maf.tsv.gz")
tumorOnly_file <- file.path(data_dir, "snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz")
expr_file <- file.path(data_dir, "gene-expression-rsem-tpm-collapsed.rds")
cn_arm_file <- file.path(data_dir, "cnv-consensus-gistic.zip")
germline_file <- file.path(input_dir, "pbta-merged-plp-variants-cbtn-plus-mb-x01.tsv")
germline_sv_file <- file.path(input_dir, "pbta_germline_svs.tsv")
tp53_alterations_file <- file.path(root_dir, "analyses",
                                  "tp53_nf1_score", "results",
                                  "tp53_altered_status.tsv")

# load histologies and subset for mb shh
hist <- read_tsv(hist_file, guess_max = 100000) %>%
  select(-molecular_subtype) %>%
  mutate(age_at_diagnosis_years = age_at_diagnosis_days/365.25)

mb_shh <- read_tsv(mb_file) %>%
  filter(grepl("SHH", molecular_subtype))

# Load gene CN data and subset for GOIs
cn_gene <- data.table::fread(cn_file) %>%
  filter(biospecimen_id %in% mb_shh$Kids_First_Biospecimen_ID) %>%
  filter(gene_symbol %in% c("MYCN", "GLI2", "CCND2", "PTEN", "PPM1D", "MDM4", "TP53"),
         status %in% c("amplification", "gain", "loss", "deep deletion"))

# load maf files and subset for GOIs
keep_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Start_Position", "End_Position",
               "HGVSp_Short", "HGVSc", "Consequence", "Variant_Classification")

mutations <- data.table::fread(maf_file) %>%
  select(all_of(keep_cols)) %>%
  filter(Hugo_Symbol %in% c("CREBBP", "DDX3X", "KMT2D", "PTCH1", "SMO", "SUFU", "TERT") |
           # for U1 mutations, keep only those at 3rd and 5th nucleotide
           grepl("RN*U1-*", Hugo_Symbol) & HGVSc %in% c("n.3A>G", "n.5A>G")) %>%
  filter(Tumor_Sample_Barcode %in% mb_shh$Kids_First_Biospecimen_ID)

tumorOnly <- data.table::fread(tumorOnly_file) %>%
  select(all_of(keep_cols)) %>%
  filter(Hugo_Symbol %in% c("CREBBP", "DDX3X", "KMT2D", "PTCH1", "SMO", "SUFU", "TERT") |
           grepl("RN*U1-*", Hugo_Symbol) & HGVSc %in% c("n.3A>G", "n.5A>G"))

mutations <- mutations %>% 
  bind_rows(tumorOnly) %>%
  filter(Tumor_Sample_Barcode %in% mb_shh$Kids_First_Biospecimen_ID)


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

# Load expr data
expr <- readRDS(expr_file)

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
  

# Create subtypes df that:
# 1) calculates average methylations score per match id
# 2) redefines SHH subtypes as alpha (3), beta (1), gamma (2), delta (4)
# 3) append molecular alterations data
# 4) attempts to classify samples without methylation classifiers into one of the four SHH subtypes using molecular and age criteria 
hist_mb_shh <- hist %>%
  filter(Kids_First_Biospecimen_ID %in% mb_shh$Kids_First_Biospecimen_ID) %>%
  select(Kids_First_Participant_ID, sample_id, match_id, age_at_diagnosis_years, dkfz_v12_methylation_subclass, dkfz_v12_methylation_subclass_score) %>%
  unique() %>%
  group_by(Kids_First_Participant_ID, sample_id, match_id, age_at_diagnosis_years) %>%
  summarize(
    dkfz_v12_methylation_subclass_collapsed = if(all(is.na(dkfz_v12_methylation_subclass))) NA_character_ 
    else paste(unique(na.omit(dkfz_v12_methylation_subclass)), collapse = "; "),
    dkfz_v12_methylation_subclass_score_mean = mean(dkfz_v12_methylation_subclass_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~replace(., is.nan(.), NA)),
         SHH_subtype = case_when(dkfz_v12_methylation_subclass_score_mean >= 0.7 & dkfz_v12_methylation_subclass_collapsed == "MB_SHH_1" ~ "SHH_beta",
                                 dkfz_v12_methylation_subclass_score_mean >= 0.7 & dkfz_v12_methylation_subclass_collapsed == "MB_SHH_2" ~ "SHH_gamma",
                                 dkfz_v12_methylation_subclass_score_mean >= 0.7 & dkfz_v12_methylation_subclass_collapsed == "MB_SHH_3" ~ "SHH_alpha",
                                 dkfz_v12_methylation_subclass_score_mean >= 0.7 & dkfz_v12_methylation_subclass_collapsed == "MB_SHH_4" ~ "SHH_delta",
                                 TRUE ~ NA_character_))

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
  # clean up subclass and scores
  mutate(dkfz_v12_methylation_subclass_score_mean = case_when(dkfz_v12_methylation_subclass_score_mean < 0.7 & 
                                                                !grepl("SHH", dkfz_v12_methylation_subclass_collapsed) ~ NA_integer_,
                                                              TRUE ~ dkfz_v12_methylation_subclass_score_mean),
         dkfz_v12_methylation_subclass_collapsed = case_when((dkfz_v12_methylation_subclass_score_mean < 0.7 | is.na(dkfz_v12_methylation_subclass_score_mean)) & 
                                                               !grepl("SHH", dkfz_v12_methylation_subclass_collapsed) ~ NA_character_,
                                                             TRUE ~ dkfz_v12_methylation_subclass_collapsed)) %>%
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

