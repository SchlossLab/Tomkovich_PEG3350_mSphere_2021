source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_histo_scores_d4.png")
b <- ggdraw() + draw_image("results/figures/5_days_PEG_histo_scores_d6.png")
c <- ggdraw() + draw_image("results/figures/5_days_PEG_tissues_PCoA.png")
d <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_tissues.png")
e <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_tissues.png")

left_panel <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol = 1, rel_widths = c(1.5, 1))
right_panel <- plot_grid(c, labels = c("C"), label_size = 12, nrow = 1)
top_panel <- plot_grid(left_panel, right_panel, nrow = 1)
bottom_panel <- plot_grid(d, e, labels = c("D", "E"), label_size = 12, nrow = 2)
plot_grid(top_panel, bottom_panel, nrow = 2, rel_widths = c(.5, 3)) +
  ggsave("results/figures/figure_3.pdf", width = 6.875, height = 9)+
  ggsave("submission/figure_3.pdf", width = 6.875, height = 9)
