library(tidyverse)
library(scales)
library(readr)
library(forcats)
library(maftools)

hist <- read_tsv("~/OpenPedCan-Project-CNH/data/histologies.tsv")
dmg_samples <- hist %>%
  filter(cancer_group %in% c("Diffuse midline glioma", "diffuse intrinsic pontine glioma") | pathology_diagnosis == "Brainstem glioma- Diffuse intrinsic pontine glioma") %>%
  filter(experimental_strategy %in% c("WGS", "WXS", "Targeted Panel"),
         tumor_descriptor %in% c("Initial CNS Tumor", "Primary Tumor"))

dmg_patient_n <- dmg_samples %>%
  select(Kids_First_Participant_ID) %>%
  unique() %>%
  nrow()

maf_to <- read_tsv("~/OpenPedCan-Project-CNH/data/snv-mutect2-tumor-only-plus-hotspots.maf.tsv.gz") %>%
  filter(Tumor_Sample_Barcode %in% dmg_samples$Kids_First_Biospecimen_ID,
         Hugo_Symbol == "ACVR1") %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Type, HGVSg, HGVSc, HGVSp, VAF, HotSpotAllele)

maf_tn <- read_tsv("~/OpenPedCan-Project-CNH/data/snv-consensus-plus-hotspots.maf.tsv.gz") %>%
  filter(Tumor_Sample_Barcode %in% dmg_samples$Kids_First_Biospecimen_ID,
         Hugo_Symbol == "ACVR1") %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Type, HGVSg, HGVSc, HGVSp, VAF, HotSpotAllele)

maf_cat <- maf_to %>%
  bind_rows(maf_tn)

maf_clin <- maf_cat %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  left_join(hist[,c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID", "age_at_event_days")]) %>%
  mutate(age_at_event_years = age_at_event_days / 365.25) %>%
  filter(!is.na(HGVSp))

# by patient tables
maf_by_patient <- maf_clin %>%
  group_by(Kids_First_Participant_ID) %>%
  summarise(
    n_biospecimens = n_distinct(Kids_First_Biospecimen_ID),
    min_age_years = min(age_at_event_years, na.rm = TRUE),
    max_age_years = max(age_at_event_years, na.rm = TRUE),
    variants = list(
      distinct(
        pick(
          Kids_First_Biospecimen_ID,
          Variant_Classification, HGVSg, HGVSc, HGVSp, HotSpotAllele,
          VAF, age_at_event_years
        )
      )
    ),
    .groups = "drop"
  )

maf_by_patient
# To view a patient's variants:
# maf_by_patient$variants[[1]]
# To expand back out:
# maf_by_patient %>% unnest(variants)

dmg_patients %>%
  filter(is.na(age_years))

# --- 1) Patient-level table for the DIPG/DMG cohort (denominator) ---
dmg_patients <- dmg_samples %>%
  select(Kids_First_Participant_ID, age_at_event_days) %>%
  mutate(age_years = ifelse(!is.na(age_at_event_days), age_at_event_days / 365.25, NA_integer_),
         age_group = case_when(age_years < 12 ~ "<12",
                               age_years >= 12 & age_years <= 18 ~ "12-18",
                               age_years > 18 ~ ">18",
                               TRUE ~ NA_character_
    ),
    age_group = factor(age_group, levels = c("<12", "12-18", ">18"))
  )

# --- 2) Patient-level mutation status (numerator) ---
acvr1_mut_patients <- maf_cat %>%
  distinct(Kids_First_Participant_ID) %>%
  mutate(ACVR1_mut = TRUE)

# --- 3) Incidence by age group ---
incidence_tbl <- dmg_patients %>%
  left_join(acvr1_mut_patients, by = "Kids_First_Participant_ID") %>%
  mutate(ACVR1_mut = coalesce(ACVR1_mut, FALSE)) %>%
  filter(!is.na(age_group)) %>%
  group_by(age_group) %>%
  summarise(
    n_patients = n(),
    n_acvr1_mut = sum(ACVR1_mut),
    incidence = n_acvr1_mut / n_patients,
    .groups = "drop"
  )

incidence_tbl

tidy_incidence <- dmg_patients %>%
  left_join(acvr1_mut_patients, by = "Kids_First_Participant_ID") %>%
  mutate(
    ACVR1_mut = coalesce(ACVR1_mut, FALSE),
    age_group = factor(age_group, levels = c("<12", "12-18", ">18"))
  ) %>%
  filter(!is.na(age_group)) %>%
  group_by(age_group) %>%
  unique() %>%
  summarise(
    total_patients = n(),
    ACVR1_mutated = sum(ACVR1_mut),
    wildtype = total_patients - ACVR1_mutated,
    incidence = ACVR1_mutated / total_patients,
    incidence_percent = percent(incidence, accuracy = 0.1),
    summary = paste0(ACVR1_mutated, "/", total_patients,
                     " (", percent(incidence, accuracy = 0.1), ")"),
    .groups = "drop"
  ) %>%
  arrange(age_group)

tidy_incidence  


## lollipops by age
# sample-level clinical annotation (maftools uses Tumor_Sample_Barcode)
clin_age <- dmg_samples %>%
  transmute(
    Tumor_Sample_Barcode = Kids_First_Biospecimen_ID,
    Kids_First_Participant_ID,
    age_years = age_at_event_days / 365.25,
    age_group = case_when(
      age_years < 12 ~ "<12",
      age_years >= 12 & age_years <= 18 ~ "12-18",
      age_years > 18 ~ ">18",
      TRUE ~ NA_character_
    )
  ) %>%
  distinct() %>%
  filter(!is.na(age_group)) %>%
  mutate(age_group = factor(age_group, levels = c("<12", "12-18", ">18")))

# make sure your maf has Tumor_Sample_Barcode (right now it's called Tumor_Sample_Barcode in maf_cat)
# and make an AA column that exists (use HGVSp stripped of "p.")
maf_for_maftools <- maf_cat %>%
  mutate(
    HGVSp_Short = str_remove(HGVSp, "^p\\.")
  )

acvr1_maf_obj <- read.maf(
  maf = maf_for_maftools,
  clinicalData = clin_age
)

pdf("~/OpenPedCan-Project-CNH/ACVR1_lollipop_by_age_group.pdf",
    width = 12, height = 5)

age_levels <- levels(clin_age$age_group)

for (ag in age_levels) {
  
  tsb_subset <- clin_age %>%
    filter(age_group == ag) %>%
    pull(Tumor_Sample_Barcode) %>%
    unique()
  
  # If there are no samples in this age bin
  if (length(tsb_subset) == 0) {
    plot.new()
    title(main = paste("ACVR1 lollipop - age group:", ag))
    text(0.5, 0.5, "No samples in this age group")
    next
  }
  
  m_sub <- tryCatch(
    subsetMaf(acvr1_maf_obj, tsb = tsb_subset),
    error = function(e) NULL
  )
  
  # If no non-syn variants (or any other subsetting error)
  if (is.null(m_sub) || nrow(m_sub@data) == 0) {
    plot.new()
    title(main = paste("ACVR1 lollipop - age group:", ag))
    text(0.5, 0.5, "No non-synonymous ACVR1 variants in this age group")
    next
  }
  
  lollipopPlot(
    maf = m_sub,
    gene = "ACVR1",
    AACol = "HGVSp_Short",
    labelPos = "all",
    repel = TRUE,
    collapsePosLabel = TRUE,
    labPosSize = 0.7,
    labPosAngle = 45,
    pointSize = 0.8,
    showMutationRate = TRUE
  )
  title(main = paste("ACVR1 lollipop - age group:", ag), line = -0.5)
}

dev.off()

