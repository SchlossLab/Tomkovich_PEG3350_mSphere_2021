source("code/utilities.R")

a <- ggdraw() + draw_image("results/figures/1_Day_PEG_Schematic.png", scale = 1.2, x = .07)
b <- ggdraw() + draw_image("results/figures/1_day_PEG_cfu.png", scale = 1.05)
c <- ggdraw() + draw_image("results/figures/1_Day_PEG_PCoA.png", scale = 1.25)
d <- ggdraw() + draw_image("results/figures/1_Day_PEG_shannon.png")
e <- ggdraw() + draw_image("results/figures/1_Day_PEG_genus_10_baselinetoD1_heatmap.png")
top_panel <- plot_grid(a, b, c, d, labels = c("A", "B", "C", "D"), label_size = 12, ncol = 2, rel_widths = c(.60, 1, 1,.60), rel_heights = c(1, 1))
bottom_panel <- plot_grid(e, labels = c("E"), label_size = 12)
plot_grid(top_panel, bottom_panel, nrow = 2, rel_heights = c(1.1, 1)) +
  ggsave("results/figures/figure_3.pdf", width = 5.25, height = 5.00)

