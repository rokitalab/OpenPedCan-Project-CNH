library(data.table)
library(tidyverse)
library(reshape2)
#library(xlsx)
library(ggrepel)
library(rprojroot)
library(ggplot2)

root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data", "v15")

if (!dir.exists("plots")) {
  dir.create("plots")
}
if (!dir.exists("results")) {
  dir.create("results")
}

histologiesfn <- file.path(data_dir, 'histologies.tsv')
#define sample id and read in only rows with that sample id
dicer_sample_id <- '7316-1746'
histologies <- fread(
  cmd = sprintf(
    "awk -F'\t' 'NR==1 || $0 ~ /(^|\\t)%s(\\t|$)/' %s",
    dicer_sample_id,
    histologiesfn
  )
)
#list of genes of interest
genes <- readLines("input/genes-of-interest.txt")

maf_file <- file.path(data_dir,"snv-consensus-plus-hotspots.maf.tsv.gz")
##this reads only the rows from the consensus maf that we want because all of the rows we want are compared to this matched normal
matched_normal <- "BS_P4Z4YP0P"
# fread with shell filtering
dicer_maf <- fread(
  cmd = sprintf(
    "zcat %s | awk 'NR==1 || $0 ~ /%s/'",
    maf_file,
    matched_normal
  )
)

#subset to genes and output select columns
dicer_maf <- dicer_maf %>%
  mutate(VAF = t_alt_count/(t_ref_count+t_alt_count)) %>%
  filter(Tumor_Sample_Barcode %in% histologies$Kids_First_Biospecimen_ID) %>% # subset to samples
  select(Hugo_Symbol, HGVSp, HGVSp_Short, VAF, Variant_Classification, Tumor_Sample_Barcode) %>%
  inner_join(histologies, by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")) %>%
  filter(Hugo_Symbol %in% genes)

##dont need this -- just for inspection 
#unique_bsaids <- unique(dicer_maf$Tumor_Sample_Barcode)

#filtered_histologies <- histologies %>%
#  filter(Kids_First_Biospecimen_ID %in% unique_bsaids)


maf_subset <- dicer_maf %>%
  filter(HGVSp_Short != "", Hugo_Symbol != "") %>%
  mutate(
    label = paste0(Hugo_Symbol, "_", HGVSp_Short),
    gene = Hugo_Symbol %in% genes
  )

maf_summary <- maf_subset %>%
  group_by(label, gene, composition, cell_line_composition, cell_line_passage) %>%
  summarize(VAF = max(VAF, na.rm = TRUE), .groups = "drop")

write.table(
  maf_summary,
  file = "results/dicer1-maf-summary.tsv",
  sep = "\t",          d
  quote = FALSE,       
  row.names = FALSE    
)
##commenting out but I may need to use later depending on feedback
#maf_wide <- maf_summary %>%
#  pivot_wider(
#    id_cols = c(label, gene),
#    names_from = c(composition, cell_line_composition, cell_line_passage),
#    values_from = VAF
#  ) %>%
#  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
#  as.data.frame()

#maf <- dicer_maf %>%
#  filter(#composition.type %in% c("Solid_Tissue", type),
#    HGVSp_Short != "",
#    Hugo_Symbol != "") %>%
#  mutate(label = paste0(Hugo_Symbol, "_", HGVSp_Short),
 #        gene = ifelse(Hugo_Symbol %in% genes, T, F)) %>%
#  group_by(label, gene, composition) %>%
#  summarize(VAF = max(VAF, na.rm = T), .groups = "drop") %>%
#  pivot_wider(id_cols = c(label, gene), names_from = composition, values_from = VAF) %>%
#  mutate_if(is.numeric, funs(replace_na(., 0))) %>%
#  as.data.frame()

# Tumor
tumor_maf <- maf_summary %>%
  filter(composition == "Solid Tissue") %>%
  select(label, VAF_tumor = VAF)

# Cell lines
cell_line_maf <- maf_summary %>%
  filter(composition == "Derived Cell Line") %>%
  select(label, VAF_cell_line = VAF, cell_line_composition, cell_line_passage)


maf_plot <- cell_line_maf %>%
  left_join(tumor_maf, by = "label") %>%
  # Replace NAs with 0 for comparison
  mutate(
    VAF_tumor = replace_na(VAF_tumor, 0),
    VAF_cell_line = replace_na(VAF_cell_line, 0),
    type = case_when(
      VAF_tumor > 0 & VAF_cell_line > 0 ~ "Common",
      VAF_tumor > 0 & VAF_cell_line == 0 ~ "Solid Tissue",
      VAF_tumor == 0 & VAF_cell_line > 0 ~ "Derived Cell Line",
      TRUE ~ NA_character_  # catches 0/0 or unexpected cases
    ))

#make sure VAFs are numeric and NA replaced by 0
maf_plot <- maf_plot %>%
  mutate(
    VAF_tumor = replace_na(VAF_tumor, 0),
    VAF_cell_line = replace_na(VAF_cell_line, 0)
  )

#Convert cell_line_passage to factor for shape mapping
maf_plot$cell_line_passage <- factor(maf_plot$cell_line_passage)

library(ggplot2)
library(ggrepel)

maf_plot$type <- factor(maf_plot$type, levels = c("Solid Tissue", "Derived Cell Line", "Common"))

# Scatter plot
p <- ggplot(maf_plot, aes(
  x = VAF_tumor, 
  y = VAF_cell_line, 
  color = cell_line_composition, 
  shape = cell_line_passage
)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Serum-based" = "firebrick3", "Serum-free" = "dodgerblue3")) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  geom_text_repel(aes(label = label), size = 3, segment.color = NA) +
  scale_shape_discrete(name = "Cell Line Passage") + 
  ggtitle("DICER-1 Patient Mutations: Tumor vs Cell Line") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold")
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("Tumor VAF") +
  ylab("Cell Line VAF") + 
  labs(color = "Cell Line Composition")


pdf(file = file.path('plots', "dicer-vaf.pdf"), width = 10, height = 8)
print(p)
dev.off()