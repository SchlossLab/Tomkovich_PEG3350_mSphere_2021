source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_stool_PCoA.png", scale = 1)
legend <- ggdraw() + draw_image("results/figures/5_days_PEG_pcoa_legend.png", scale = 1)
legend<- plot_grid(legend, NULL, NULL, labels = NULL, label_size = 12, nrow = 1, rel_heights = c(1.5, 1), rel_widths = c(1, 0.3))
b <- ggdraw() + draw_image("results/figures/5_days_PEG_bc_combined_baseline.png", scale = 1)
c <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_stools_median.png", scale = 1)
d <- ggdraw() + draw_image("results/figures/5_days_PEG_genera_impacted_by_PEG_barbell_v2.png")
e <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_stools.png", scale =0.95)

top_panel <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, nrow = 1, rel_widths = c(1, 1.2))
middle_panel <- plot_grid(c, labels = c("C"), label_size = 12, ncol = 1)
b_right_panel <- plot_grid(e, NULL, labels = c("E", ""), label_size = 12, ncol = 1)
bottom_panel <- plot_grid(d, b_right_panel, labels = c("D", ""), label_size = 12, nrow = 1, rel_heights = c(1.2, 1))
plot_grid(top_panel, legend, middle_panel, bottom_panel, ncol = 1, rel_heights = c(1, 0.2, 1, 2)) +
  ggsave("results/figures/figure_2.pdf", width = 6.875, height = 9)+
  ggsave("submission/figure_2.pdf", width = 6.875, height = 9)

