source("code/utilities.R")

a <- ggdraw() + draw_image("results/pictures/1_Day_PEG_Schematic.png", scale = 1.5, x = .06, y = -.1)
b <- ggdraw() + draw_image("results/figures/1_day_PEG_cfu.png")
c <- ggdraw() + draw_image("results/figures/1_Day_PEG_PCoA.png", x = -.1, y = -.05)
d <- ggdraw() + draw_image("results/figures/1_Day_PEG_shannon.png")
e <- ggdraw() + draw_image("results/figures/1_Day_PEG_genus_10_PTtoD1_heatmap.png")
left_panel <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol = 1)
right_panel <- plot_grid(c, d, labels = c("C", "D"), label_size = 12, ncol = 1)
bottom_panel <- plot_grid(e, labels = c("E"), label_size = 12, ncol = 1)
plot_grid(left_panel, right_panel, bottom_panel) +
  ggsave("results/figures/figure_3.pdf", width=6.875, height=4.73)

