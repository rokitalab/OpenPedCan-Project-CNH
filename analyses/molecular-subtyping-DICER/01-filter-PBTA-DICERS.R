library(tidyverse)
library(dplyr)
library(readr)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

histologies_df <- read_tsv(file.path(data_dir, "histologies-base.tsv"), 
                           guess_max = 100000)


results_dir <- "results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

biospecimen_id_fn <- file.path(results_dir, "biospecimen_ids_PBTA_no_methylation.tsv")

##want to filter out kids first pids any that have methylation classifcations above .8 that are not DICER

histologies_filtered <- histologies_df %>%
  filter(
    cohort == "PBTA",                                        # keep PBTA only
    !(dkfz_v12_methylation_subclass != "CNS_SARC_DICER" & 
        dkfz_v12_methylation_subclass_score >= 0.8)            # remove flagged patients
  )



pbta_biospecimen_ids <- histologies_filtered %>%
  pull(Kids_First_Biospecimen_ID)

##this gives use ids to use for reading the SNV calls. 
write.table(pbta_biospecimen_ids,
            file = biospecimen_id_fn,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
#incase I need to use later
write_tsv(histologies_filtered, "results/histologies_filtered.tsv")

##read snvs in chunks 
chunks <- list()  # list to store each chunk

callback <- function(x, pos) {
  chunks[[length(chunks) + 1]] <<- x %>%
    filter(Hugo_Symbol == "DICER1")
}

read_tsv_chunked(
  file.path(data_dir, "snv-consensus-plus-hotspots.maf.tsv.gz"),
  callback = DataFrameCallback$new(callback),
  chunk_size = 50000
)

# combine all chunks into one tibble
dicer1_mut <- bind_rows(chunks)



filter_dicer1_variants_types <- function(dicer1_mut_df) {
  dicer1_mut_df %>%
    filter(
      Tumor_Sample_Barcode %in% pbta_biospecimen_ids,
      Variant_Classification %in% c(
        "Missense_Mutation",
        "Nonsense_Mutation",
        "Frame_Shift_Ins",
        "Frame_Shift_Del"
      )
    )
}



filter_dicer1_variants_types <- function(dicer1_mut_df, pbta_biospecimen_ids) {
  dicer1_mut_df %>%
      filter(
        Tumor_Sample_Barcode %in% pbta_biospecimen_ids,
        Variant_Classification %in% c(
          "Missense_Mutation",
          "Nonsense_Mutation",
          "Frame_Shift_Ins",
          "Frame_Shift_Del"
        )
      )
}


dicer1_filtered <- filter_dicer1_variants_types(dicer1_mut, pbta_biospecimen_ids)

write_tsv(dicer1_filtered, "results/dicer1_variants_filtered_final.tsv")


