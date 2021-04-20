source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_tissues.png")
b <- ggdraw() + draw_image("results/figures/5_days_PEG_tissues_PCoA.png")
c <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_tissues.png")
middle_panel <- plot_grid(b, c, labels = c("B", "C"), label_size = 12, nrow = 1, rel_widths = c(.6, 1))
d <- ggdraw() + draw_image("results/figures/5_days_PEG_histo_scores_d4.png")
e <- ggdraw() + draw_image("results/figures/5_days_PEG_histo_scores_d6.png")
bottom_panel <- plot_grid(d, e, labels = c("D", "E"), label_size = 12, nrow = 1, rel_heights = c(1.2, 1))

plot_grid(a, middle_panel, bottom_panel, labels = c("A", "", ""), nrow = 3) +
  ggsave("results/figures/figure_3.pdf", width = 6.875, height = 9)+
  ggsave("submission/figure_3.pdf", width = 6.875, height = 9)
