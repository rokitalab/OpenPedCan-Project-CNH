# Plot theme for the OpenPedCan figures: the look of the manuscript theme
# (figures/manuscript_OPC/utils/theme_for_plot.R) on a white background, built
# on ggplot2's theme_bw() so it needs no extra package.
#
# Usage:
#   source(file.path(root_dir, "figures", "scripts", "figure_theme.R"))
#   ggplot(...) + theme_Publication()
#
# Panels, the plot, legend keys and facet labels are white, the axes are black
# lines, there is no panel border, and the major grid lines are light grey
# (use `+ theme(panel.grid.major = element_blank())` to drop them).

theme_Publication <- function(base_size = 15, base_family = "Helvetica") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2),
                                         hjust = 0.5),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#f0f0f0"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
      axis.title.x = ggplot2::element_text(vjust = -0.2),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.ticks = ggplot2::element_line(),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.direction = "vertical",
      legend.key.size = ggplot2::unit(0.5, "cm"),
      legend.margin = ggplot2::margin(5, 5, 5, 5),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::unit(c(10, 5, 5, 10), "mm"),
      strip.background = ggplot2::element_rect(colour = "#f0f0f0", fill = "white"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
