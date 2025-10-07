# Stephanie J. Spielman and Jaclyn Taroni for ALSF CCDL 2020
#
# This script subsets the files required for subtyping non-MB and non-ATRT
# embryonal tumors. The samples that were subset in 
# [`01-samples-to-subset.Rmd`](./01-samples-to-subset.Rmd), based on specific
# conditions outlined in that notebook, will be included here.

library(tidyverse)
library(readr)
library(dplyr)

#### Directories ---------------------------------------------------------------

# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "analyses", "molecular-subtyping-embryonal")

# Set path to directory where we will save the subset files
subset_dir <- file.path(module_dir, "subset-files")
if (!dir.exists(subset_dir)) {
  dir.create(subset_dir)
}

# Results directory
results_dir <- file.path(module_dir, "results")

#### Read in files -------------------------------------------------------------

# File that contains the relevant biospecimen identifiers -- this was generated
# in 01-samples-to-subset
subset_id <-
  read_tsv(file.path(results_dir,
                     "biospecimen_ids_embryonal_subtyping.tsv")) %>%
  pull(Kids_First_Biospecimen_ID)

# There are relevant samples in RNA-seq data
rnaseq_collapsed <-
  read_rds(file.path(data_dir,
                     "gene-expression-rsem-tpm-collapsed.rds"))

# structural variant for BCOR tandem duplications
manta_sv_df <- data.table::fread(file.path(data_dir, "sv-manta.tsv.gz"))

#### Filter and process expression data ----------------------------------------

# In this particular case, we're going to z-score before subsetting because
# there are few samples -- so the values that end up in the final table will
# be in the context of the entire cohort

filter_process_expression <- function(expression_mat) {
  # This function takes an expression matrix where the columns are samples
  # and the rows are genes, where the gene identifiers are rownames.
  # It log transforms the matrix, z-scores the genes, and then filters the
  # samples to only those that are in `subset_id`.
  # It returns matrix where the genes are columns and the rows are samples.
  # ONLY INTENDED TO BE USED IN THIS CONTEXT!

  # log2(x + 1) transform the expression matrix
  log_expression <- log2(expression_mat + 1)
  # Scale the gene values -- scale() works on the columns, hence the transpose
  z_scored_expression <- scale(t(log_expression),
                               center = TRUE,
                               scale = TRUE)

  # return z-scored expression data
  return(z_scored_expression)
}

# filter for embryonal samples and then z-score
rnaseq_collapsed <- rnaseq_collapsed %>%
  dplyr::select(intersect(subset_id, colnames(rnaseq_collapsed)))
filter_process_expression(rnaseq_collapsed) %>%
  write_rds(file.path(subset_dir, "embryonal_zscored_exp.rds"))

#### Structural variant data ---------------------------------------------------

manta_sv_df %>%
  filter(Kids.First.Biospecimen.ID.Tumor %in% subset_id) %>%
  write_tsv(file.path(subset_dir, "embryonal_manta_sv.tsv"))

rm(manta_sv_df)

#### SNV analysis

 




# Storage list for collected chunks
chunks <- list()

# read in dicer 1 mutants that have tumor sample barcode in subset id
callback <- function(x, pos) {
  chunks[[length(chunks) + 1]] <<- x %>%
    filter(
      Hugo_Symbol == "DICER1",
      Tumor_Sample_Barcode %in% subset_id,

    )
}


read_tsv_chunked(
  file.path(data_dir, "snv-consensus-plus-hotspots.maf.tsv.gz"),
  callback = DataFrameCallback$new(callback),
  chunk_size = 50000
)

# combine all chunks into one tibble
dicer1_mut <- bind_rows(chunks)
#filter variants to only deliterious, and write only the biospecimen id and the HGVsc to output file 
dicer1_mut %>%
  filter(
    Variant_Classification %in% c(
      "Missense_Mutation",
      "Nonsense_Mutation",
      "Frame_Shift_Ins",
      "Frame_Shift_Del"
    )
  ) %>%
  select(
    `Kids_First_Biospecimen_ID` = Tumor_Sample_Barcode,
    `DICER-1_Mutation` = HGVSc
  ) %>%
  write_tsv(file.path(subset_dir, "dicer1_mutations.tsv"))

