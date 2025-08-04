# ============================================================================================================
# Script Name: immune_cell_medulloblastoma_analysis.R
# Purpose: Generate immune cell comparison between medulloblastoma subtypes
# Author: Daniel Emerson
# Date: 2025-08-02
#
# Summary:
# Group3 is most significantly different from WNT, group2, and SHH for immune fraction.
# (smallest pvalue comparison (both Mann Whitney and KS test on both total immune fraction and cell
# subtypes (e.g. monocytes). Greatest shift in median of distribution).
#
# I tested both total (sum of celltypes and all separately) and looked at particular celltype such as monocytes.
# I generated distributions of all 6 combinations, performed non parametric stats tests and volcano plots 
#  comparing celltype fraction enrichment for group3 vs every other type. I used quantiseq, which provides 
# absolute fraction as opposed to an arbitrary enrichment score.  Caveat: not able to characterize as many celltypes 
# which are lost as "uncharacterized" with less sensitivity, which may impact analysis.
# 
# Output files:
# distributions_immune_cell_subtype_medulloblastoma.png - pairwise distribution comparison total immune fraction
# distributions_monocyte_cell_subtype_medulloblastoma.png - pairwise distribution comparison monocyte fraction
# heatmap_immune_cell_subtype_medulloblastoma.png - heatmap of mann whitney u and ks pairwise for total immune 
# volcano_group3_vs_rest_immune_cell_type.png - pval vs logFC for group3 vs all three (rest) for all celltypes
# volcano_group3_vs_group4_immune_cell_type.png - pval vs logFC for group3 vs group4 for all celltypes
# volcano_group3_vs_SHH_immune_cell_type.png - pval vs logFC for group3 vs SHH for all celltypes
# volcano_group3_vs_WNT_immune_cell_type.png - pval vs logFC for group3 vs WNT for all celltypes
# ============================================================================================================


# ---- Load Libraries ----

library(ggplot2)
library(dplyr)
library(pheatmap)
library(reshape2)
library(cowplot)

setwd("./")


# ---- load and filter tumor file ----
#extract xcell based deconvolution of celltypes
#deconv_xcell_df <- readRDS("results/xcell_output.rds")

#extract quantiseq based deconvolution of celltypes
deconv_quantiseq_output_df <- readRDS("results/quantiseq_output.rds")
deconv_quantiseq_output_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$cancer_group == "Medulloblastoma",]

#remove uncharacterized celltype fraction
deconv_quantiseq_output_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$cell_type != "uncharacterized cell",]

#create total fraction sum of immune cells per subtype per individual
total_frac_immune_df <-  deconv_quantiseq_output_df  %>% group_by(Kids_First_Biospecimen_ID) %>%
  reframe(fraction = sum(fraction), molecular_subtype = molecular_subtype,cancer_group = cancer_group)
total_frac_immune_df <- unique(total_frac_immune_df)                             


#filter for Medulloblastoma tumor and subtype

#went with total fraction sum since massive zero fraction for cell types distorts outcome of distribution
# and corresponding medians.  Medians became close to zero as a result.

MB_group4_df <- total_frac_immune_df[total_frac_immune_df$cancer_group == "Medulloblastoma" & 
                                       total_frac_immune_df$molecular_subtype == "MB, Group4",]

#alternative test with each type seprate in total distribution
#MB_group4_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$cancer_group == "Medulloblastoma" &
#                                             deconv_quantiseq_output_df$molecular_subtype == "MB, Group4",]

MB_group3_df <- total_frac_immune_df[total_frac_immune_df$cancer_group == "Medulloblastoma" & 
                                      total_frac_immune_df$molecular_subtype == "MB, Group3",]

#alternative test with each type seprate in total distribution
#MB_group3_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$cancer_group == "Medulloblastoma" &
#                                             deconv_quantiseq_output_df$molecular_subtype == "MB, Group3",]

MB_SHH_df <- total_frac_immune_df[total_frac_immune_df$cancer_group == "Medulloblastoma" &      
                                    total_frac_immune_df$molecular_subtype == "MB, SHH",]

#alternative test with each type seprate in total distribution
#MB_SHH_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$cancer_group == "Medulloblastoma" & 
#                                          deconv_quantiseq_output_df$molecular_subtype == "MB, SHH",]

MB_WNT_df <- total_frac_immune_df[total_frac_immune_df$cancer_group == "Medulloblastoma" & 
                                    total_frac_immune_df$molecular_subtype == "MB, WNT",]

#alternative test with each type seprate in total distribution
#MB_WNT_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$cancer_group == "Medulloblastoma" & 
#                                          deconv_quantiseq_output_df$molecular_subtype == "MB, WNT",]


# ---- pairwise statistical tests for total immune cell fraction by type  ----
# perform pairwise test (6 total) comparing total immune cell fractions distribution for 4 conditions
# (WNT, group4,group3,SHH) using non-parametric tests Mann-Whitney U and KS
# I performed KS and Mann Whitney U as non parametric rank based tests (no assumptions about
# distribution type/parameters).  I looked for consistancy in trends. 

WNT_v_group4 <- wilcox.test(MB_WNT_df$fraction,MB_group4_df$fraction,alternative = "two.sided")
WNT_v_group4_ks <- ks.test(MB_WNT_df$fraction,MB_group4_df$fraction, alternative = "two.sided")

WNT_v_group3 <- wilcox.test(MB_WNT_df$fraction,MB_group3_df$fraction,alternative = "two.sided")
WNT_v_group3_ks <- ks.test(MB_WNT_df$fraction,MB_group3_df$fraction, alternative = "two.sided")

WNT_v_SHH <- wilcox.test(MB_WNT_df$fraction,MB_SHH_df$fraction,alternative = "two.sided")
WNT_v_SHH_ks <- ks.test(MB_WNT_df$fraction,MB_SHH_df$fraction, alternative = "two.sided")

group3_v_group4 <- wilcox.test(MB_group3_df$fraction,MB_group4_df$fraction,alternative = "two.sided")
group3_v_group4_ks <- ks.test(MB_group3_df$fraction,MB_group4_df$fraction, alternative = "two.sided")

group3_v_SHH <- wilcox.test(MB_group3_df$fraction,MB_SHH_df$fraction,alternative = "two.sided")
group3_v_SHH_ks <- ks.test(MB_group3_df$fraction,MB_SHH_df$fraction, alternative = "two.sided")

group4_v_SHH <- wilcox.test(MB_group4_df$fraction,MB_SHH_df$fraction,alternative = "two.sided")
group4_v_SHH_ks <- ks.test(MB_group4_df$fraction,MB_SHH_df$fraction, alternative = "two.sided")

# ---- table of pairwise pvalues for total immune cell fraction ----

sig_df <- data.frame(
  mwu = c(WNT_v_group4$p.value, WNT_v_group3$p.value, WNT_v_SHH$p.value, 
          group3_v_group4$p.value, group3_v_SHH$p.value,group4_v_SHH$p.value),
  ks = c(WNT_v_group4_ks$p.value, WNT_v_group3_ks$p.value, WNT_v_SHH_ks$p.value, 
         group3_v_group4_ks$p.value, group3_v_SHH_ks$p.value,group4_v_SHH_ks$p.value)
)

rownames(sig_df) <- c("WNT_v_group4", "WNT_v_group3", "WNT_v_SHH","group3_v_group4",
                      "group3_v_SHH","group4_v_SHH")

#It appears that group3 distribution is different from all the other groups
#Pairwise p-values greatest with group3 versus every other group

heatmap <- pheatmap(as.matrix(-log10(sig_df)), 
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           display_numbers = TRUE,
           color = colorRampPalette(c("red", "white", "blue"))(50),
           main = "-log10(P-value) Heatmap")

ggsave("plots/heatmap_immune_cell_subtype_medulloblastoma.png", plot = heatmap,width = 8, height = 6, units = "in", dpi = 300)

# ---- plot total immune cell distributions per subtype ---- 

comparisons <- new.env()
comparisons[["WNT_v_group4"]] <- c("MB, WNT","MB, Group4")
comparisons[["WNT_v_group3"]] =  c("MB, WNT","MB, Group3")
comparisons[["WNT_v_SHH"]] = c("MB, WNT","MB, SHH")
comparisons[["group3_v_group4"]] = c("MB, Group3","MB, Group4")
comparisons[["group3_v_SHH"]] = c("MB, Group3","MB, SHH")
comparisons[["group4_v_SHH"]] = c("MB, Group4","MB, SHH")


plots <- list()
color_dict <- c(
  "MB, WNT" = "skyblue",
  "MB, Group4" = "purple",
  "MB, Group3" = "salmon",
  "MB, SHH" = "green"
)  
for (comparison in names(comparisons)){

  comparison_immune_df <- total_frac_immune_df %>% filter(molecular_subtype %in% comparisons[[comparison]])
  #comparison_immune_df <- deconv_quantiseq_output_df %>% filter(molecular_subtype %in% comparisons[[comparison]])
  
  #adding medians for each comparison to plots
  group_medians <- comparison_immune_df %>%
    group_by(molecular_subtype) %>%
    summarise(median_val = median(fraction))
  
  molecular_subtype_instances <- unique(comparison_immune_df$molecular_subtype)
  
  plots[[length(plots) + 1]]  <- ggplot(comparison_immune_df, 
    aes(x = fraction, fill = molecular_subtype)) +
    geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.005) +
    geom_vline(data = group_medians, 
      aes(xintercept = median_val),
      linetype = "solid", size = 1) +
    scale_fill_manual(values = color_dict[molecular_subtype_instances]) +
    geom_text(x = group_medians$median_val[1]+0.02, y = 5, label = paste0("median: ",signif(group_medians$median_val[1],4)),
              color = color_dict[group_medians$molecular_subtype[1]], size = 2) +
    geom_text(x = group_medians$median_val[2]+0.02, y = 4, label = paste0("median: ",signif(group_medians$median_val[2],4)),
              color = color_dict[group_medians$molecular_subtype[2]], size = 2) +
    geom_text(x = 0.5, y = 5, label = paste0("p-value: ",signif(sig_df[rownames(sig_df) == comparison,][1],6)),
              color = "black", size = 3) +
    xlim(0, 0.75) +
    ylim(0, 6) +
    labs(x = "total immune cell fraction", y = "frequency") +
    theme_minimal() + 
    theme(
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "white", color = NA)
      )

}

immune_distributions <- do.call(plot_grid, c(plots, ncol = 1, align = "v"))
ggsave("plots/distributions_immune_cell_subtype_medulloblastoma.png", plot = immune_distributions,width = 8, height = 6, units = "in", dpi = 300)

#Most significant difference is with group 3


# ---- volcano plot of immune cell types in group 3 vs rest ----

pval_cell_types <- c()
fc_cell_types <- c()
cell_types <- c()
for (cell_type in unique(deconv_quantiseq_output_df$cell_type)){
  cell_type_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$cell_type == cell_type,]
  group3_cell_type_df <- cell_type_df[cell_type_df$molecular_subtype == "MB, Group3",]
  other_cell_type_df <-  cell_type_df[cell_type_df$molecular_subtype != "MB, Group3",]
  #other_cell_type_df <-  cell_type_df[cell_type_df$molecular_subtype == "MB, WNT",]
  #other_cell_type_df <-  cell_type_df[cell_type_df$molecular_subtype == "MB, Group4",]
  #other_cell_type_df <-  cell_type_df[cell_type_df$molecular_subtype == "MB, SHH",]
  
  group3_v_other_wmu <- wilcox.test(group3_cell_type_df$fraction,other_cell_type_df$fraction,alternative = "two.sided")
  pval_cell_types  <- c(pval_cell_types,-log10(group3_v_other_wmu$p.value))
  group3_v_other_fc <- mean(group3_cell_type_df$fraction)/mean(other_cell_type_df$fraction)
  fc_cell_types <- c(fc_cell_types,log2(group3_v_other_fc))
  cell_types <- c(cell_types,cell_type)
}

volcano_df <- data.frame(
  cell_type = cell_types,
  log2FC = fc_cell_types,
  pvalue = pval_cell_types
)


# Volcano plot
volcplot <- ggplot(volcano_df, aes(x = log2FC, y = pvalue)) +
  geom_point(aes(color = 'red')) +
  geom_text(data = volcano_df,aes(label = cell_type),size = 4) +
  #xlim(-5, 5) +
  #labs(title = "group 3 immune cell fraction change vs rest",x = "log2(mean(group3_fraction)/mean(rest_fraction))", y = "-log10(p-value)") +
  #labs(title = "group 3 immune cell fraction change vs WNT",x = "log2(mean(group3_fraction)/mean(rest_fraction))", y = "-log10(p-value)") +
  #labs(title = "group 3 immune cell fraction change vs group4",x = "log2(mean(group3_fraction)/mean(rest_fraction))", y = "-log10(p-value)") +
  labs(title = "group 3 immune cell fraction change vs SHH",x = "log2(mean(group3_fraction)/mean(rest_fraction))", y = "-log10(p-value)") +
  theme(legend.position = "none")

ggsave("plots/volcano_group3_vs_rest_immune_cell_type.png", plot = volcplot,width = 8, height = 6, units = "in", dpi = 300)
#ggsave("plots/volcano_group3_vs_WNT_immune_cell_type.png", plot = volcplot,width = 8, height = 6, units = "in", dpi = 300)
#ggsave("plots/volcano_group3_vs_group4_immune_cell_type.png", plot = volcplot,width = 8, height = 6, units = "in", dpi = 300)
#ggsave("plots/volcano_group3_vs_SHH_immune_cell_type.png", plot = volcplot,width = 8, height = 6, units = "in", dpi = 300)

#There appears to be an upregulation in Tcells (CD8 and CD4) and Monocytes and downregulation in rest of cells
#based on fraction of celltype found


# ---- extra pairwise statistical test for Monocyte change ----
# A valid criticism of above total sum  approach is that certain immune celltypes may upregulate and downregulate.
# Thus, there maybe a cancellation effect not caught by total sum. I chose a test focused on celltype monocyt.
# This gave largely same trends.

MB_group3_monoctye_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$molecular_subtype == 'MB, Group3' &
                             deconv_quantiseq_output_df$cell_type == "Monocyte",]       
MB_group4_monoctye_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$molecular_subtype == 'MB, Group4' &
                             deconv_quantiseq_output_df$cell_type == "Monocyte",]   
MB_WNT_monoctye_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$molecular_subtype == 'MB, WNT' &
                                                      deconv_quantiseq_output_df$cell_type == "Monocyte",]
MB_SHH_monoctye_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$molecular_subtype == 'MB, SHH' &
                                                      deconv_quantiseq_output_df$cell_type == "Monocyte",]

WNT_v_group4_monocyte <- wilcox.test(MB_WNT_monoctye_df$fraction,MB_group4_monoctye_df$fraction,
                                     alternative = "two.sided")
WNT_v_group3_monocyte <- wilcox.test(MB_WNT_monoctye_df$fraction,MB_group3_monoctye_df$fraction,
                                     alternative = "two.sided")
WNT_v_SHH_monocyte <- wilcox.test(MB_WNT_monoctye_df$fraction,MB_SHH_monoctye_df$fraction,
                                     alternative = "two.sided")
group3_v_group4_monocyte <- wilcox.test(MB_group3_monoctye_df$fraction,MB_group4_monoctye_df$fraction,
                                  alternative = "two.sided")
group3_v_SHH_monocyte <- wilcox.test(MB_group3_monoctye_df$fraction,MB_SHH_monoctye_df$fraction,
                                        alternative = "two.sided")
group4_v_SHH_monocyte <- wilcox.test(MB_group4_monoctye_df$fraction,MB_SHH_monoctye_df$fraction,
                                     alternative = "two.sided")

sig_df <- data.frame(
  mwu = c(WNT_v_group4_monocyte$p.value, WNT_v_group3_monocyte$p.value, WNT_v_SHH_monocyte$p.value, 
          group3_v_group4_monocyte$p.value, group3_v_SHH_monocyte$p.value,group4_v_SHH_monocyte$p.value)
)

rownames(sig_df) <- c("WNT_v_group4", "WNT_v_group3", "WNT_v_SHH","group3_v_group4",
                      "group3_v_SHH","group4_v_SHH")

# ---- plot monocyte distributions per subtype ----

comparisons <- new.env()
comparisons[["WNT_v_group4"]] <- c("MB, WNT","MB, Group4")
comparisons[["WNT_v_group3"]] =  c("MB, WNT","MB, Group3")
comparisons[["WNT_v_SHH"]] = c("MB, WNT","MB, SHH")
comparisons[["group3_v_group4"]] = c("MB, Group3","MB, Group4")
comparisons[["group3_v_SHH"]] = c("MB, Group3","MB, SHH")
comparisons[["group4_v_SHH"]] = c("MB, Group4","MB, SHH")


plots <- list()
for (comparison in names(comparisons)){
  
  deconv_quantiseq_output_monocyte_df <- deconv_quantiseq_output_df[deconv_quantiseq_output_df$cell_type == "Monocyte",]
  comparison_immune_df <- deconv_quantiseq_output_monocyte_df %>% 
    filter(molecular_subtype %in% comparisons[[comparison]])
  
  #adding medians for each comparison to plots
  group_medians <- comparison_immune_df %>%
    group_by(molecular_subtype) %>%
    summarise(median_val = median(fraction))
  
  molecular_subtype_instances <- unique(comparison_immune_df$molecular_subtype)
  
  plots[[length(plots) + 1]]  <- ggplot(comparison_immune_df, 
                                        aes(x = fraction, fill = molecular_subtype)) +
    geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.005) +
    geom_vline(data = group_medians, 
               aes(xintercept = median_val),
               linetype = "solid", size = 1) +
    scale_fill_manual(values = color_dict[molecular_subtype_instances]) +
    geom_text(x = group_medians$median_val[1]+0.02, y = 5, label = paste0("median: ",signif(group_medians$median_val[1],4)),
              color = color_dict[group_medians$molecular_subtype[1]], size = 2) +
    geom_text(x = group_medians$median_val[2]+0.02, y = 4, label = paste0("median: ",signif(group_medians$median_val[2],4)),
              color = color_dict[group_medians$molecular_subtype[2]], size = 2) +
    geom_text(x = 0.5, y = 5, label = paste0("p-value: ",signif(sig_df[rownames(sig_df) == comparison,][1],6)),
              color = "black", size = 3) +
    xlim(0, 0.75) +
    ylim(0, 6) +
    labs(x = "monocyte cell fraction", y = "frequency") +
    theme_minimal() + 
    theme(
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
}

immune_distributions <- do.call(plot_grid, c(plots, ncol = 1, align = "v"))
ggsave("plots/distributions_monocyte_cell_subtype_medulloblastoma.png", plot = immune_distributions,width = 8, height = 6, units = "in", dpi = 300)



