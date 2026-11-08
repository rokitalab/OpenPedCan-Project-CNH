################################################################################
# plot-palette.R
#
# The plot groups of figures/palettes/plot-mapping.tsv next to their colors (with
# the number of tumors of the integrated histologies table), and the closest pair
# of colors (CIELAB distance), to see the palette.
#
# Usage: Rscript --vanilla figures/scripts/plot-palette.R
# Output: figures/plots/plot-group-palette.pdf
################################################################################

suppressMessages(library(tidyverse))
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)
m <- read_tsv("figures/palettes/plot-mapping.tsv", show_col_types=FALSE, na=character(), col_types=cols(.default="c"))
groups <- m |> count(plot_group, plot_group_hex) |> arrange(desc(n)) |> distinct(plot_group, .keep_all=TRUE)
dups <- m |> distinct(plot_group, plot_group_hex) |> count(plot_group) |> filter(n>1)
cat("plot groups:", nrow(groups), "; groups with more than one hex:", nrow(dups), "\n")
if (nrow(dups)>0) print(dups)
h <- read_tsv("analyses/molecular-subtyping-integrate/results/histologies.tsv", show_col_types=FALSE, guess_max=1e5) |> filter(sample_type=="Tumor") |>
  mutate(tk=coalesce(match_id,Kids_First_Biospecimen_ID), cancer_group=replace_na(cancer_group,"NA")) |> distinct(tk,cohort,broad_histology,cancer_group)
tum <- h |> inner_join(m |> select(cohort,broad_histology,cancer_group,plot_group), by=c("cohort","broad_histology","cancer_group")) |> count(plot_group, name="tumors")
groups <- groups |> left_join(tum, by="plot_group") |> mutate(tumors=replace_na(tumors,0L)) |> arrange(plot_group)
# perceptual distance between hexes (CIELAB)
lab <- convertColor(t(col2rgb(groups$plot_group_hex))/255, "sRGB", "Lab")
d <- as.matrix(dist(lab)); diag(d) <- Inf
nn <- tibble(plot_group=groups$plot_group, nearest=groups$plot_group[apply(d,1,which.min)], deltaE=apply(d,1,min)) |> arrange(deltaE)
print(head(nn, 10)); 
lum <- function(hex){ r <- col2rgb(hex)/255; r <- ifelse(r<=0.03928, r/12.92, ((r+0.055)/1.055)^2.4); 0.2126*r[1,]+0.7152*r[2,]+0.0722*r[3,] }
groups <- groups |> mutate(ink = ifelse(lum(plot_group_hex) > 0.35, "#222222", "#ffffff"),
                           plot_group = factor(plot_group, levels = rev(plot_group)))
p <- ggplot(groups, aes(y = plot_group)) +
  geom_tile(aes(x = 0.5, fill = plot_group_hex), width = 1, height = 0.82, colour = "#d0d0d0", linewidth = 0.3) +
  geom_text(aes(x = 0.5, label = plot_group_hex, colour = ink), size = 3.2, family = "mono") +
  geom_text(aes(x = 1.12, label = plot_group), hjust = 0, size = 3.6, colour = "#222222") +
  geom_text(aes(x = 3.6, label = paste0(format(tumors, big.mark = ","), " tumors")), hjust = 1, size = 3.2, colour = "#666666") +
  scale_fill_identity() + scale_colour_identity() +
  scale_x_continuous(limits = c(0, 3.7), expand = c(0, 0)) +
  labs(title = "Plot group palette", subtitle = paste0(nrow(groups), " plot groups of figures/palettes/plot-mapping.tsv (tumor counts from the integrated histologies)"), x = NULL, y = NULL) +
  theme_void(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(colour = "#555555", size = 9), plot.margin = margin(12, 16, 12, 12), plot.background = element_rect(fill = "white", colour = NA))
ggsave("figures/plots/plot-group-palette.pdf", p, width = 7.5, height = 0.32 * nrow(groups) + 1.4, limitsize = FALSE)
