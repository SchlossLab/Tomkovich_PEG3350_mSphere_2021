source("code/utilities.R")

a <- ggdraw() + draw_image("results/pictures/1_Day_PEG_Schematic.png", scale = 1, x = .07)
b <- ggdraw() + draw_image("results/figures/1_day_PEG_cfu.png")
c <- ggdraw() + draw_image("results/figures/1_Day_PEG_PCoA.png", scale = 1.10)
d <- ggdraw() + draw_image("results/figures/1_Day_PEG_shannon.png")
e <- ggdraw() + draw_image("results/figures/1_Day_PEG_genus_10_baselinetoD1_heatmap.png")
top_panel <- plot_grid(a, b, c, d, labels = c("A", "B", "C", "D"), label_size = 12, ncol = 2, rel_widths = c(.55, 1, 1,.55), rel_heights = c(1, 1.1))
bottom_panel <- plot_grid(e, labels = c("E"), label_size = 12)
plot_grid(top_panel, bottom_panel, nrow = 2, rel_heights = c(1, .75)) +
  ggsave("results/figures/figure_3.pdf", width=6.875, height=4.73)

