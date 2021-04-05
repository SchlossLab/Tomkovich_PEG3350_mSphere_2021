source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_stool_PCoA.png", scale = 1)
b <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_stools.png", scale = 1.1)
c <- ggdraw() + draw_image("results/figures/5_days_PEG_genera_impacted_by_PEG.png")
d <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_stools.png", scale =0.95)

left_panel <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol = 1, rel_widths = c(1.5, 1))
right_panel <- plot_grid(c, labels = c("C"), label_size = 12, nrow = 1)
top_panel <- plot_grid(left_panel, right_panel, nrow = 1)
bottom_panel <- plot_grid(d, labels = c("D"), label_size = 12, ncol = 1)
plot_grid(top_panel, bottom_panel, nrow = 2) +
  ggsave("results/figures/figure_2.pdf", width = 6.875, height = 9)+
  ggsave("submission/figure_2.pdf", width = 6.875, height = 9)

