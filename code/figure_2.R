source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_stool_PCoA.png", scale = 1)
b <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_stools.png", scale = 1)
c <- ggdraw() + draw_image("results/figures/5_days_PEG_genera_impacted_by_PEG.png")
d <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_stools.png", scale =0.95)

top_panel <- plot_grid(a, NULL, labels = c("A", ""), label_size = 12, nrow = 1, rel_widths = c(1, 1.5))
middle_panel <- plot_grid(b, labels = c("B"), label_size = 12, ncol = 1)
b_right_panel <- plot_grid(d, NULL, labels = c("D", ""), label_size = 12, ncol = 1)
bottom_panel <- plot_grid(c, b_right_panel, labels = c("C", ""), label_size = 12, nrow = 1, rel_heights = c(1.2, 1))
plot_grid(top_panel, middle_panel, bottom_panel, ncol = 1, rel_heights = c(1, 1, 2)) +
  ggsave("results/figures/figure_2.pdf", width = 6.875, height = 9)+
  ggsave("submission/figure_2.pdf", width = 6.875, height = 9)

