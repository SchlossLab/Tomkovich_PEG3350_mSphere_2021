source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_stool_PCoA_mock.png")
a_legend <- ggdraw() + draw_image("results/figures/5_days_PEG_pcoa_mock_legend.png")
a <- plot_grid(a, a_legend, ncol = 1, rel_heights = c(1, .2))
d <- ggdraw() + draw_image("results/figures/5_days_PEG_tissue_PCoA_mock.png")
d_legend <- ggdraw() + draw_image("results/figures/5_days_PEG_pcoa_mock_alpha_legend_tissues.png")
d <- plot_grid(d, d_legend, ncol = 1, rel_heights = c(1, .1))
top_panel <- plot_grid(a, d, labels = c("A", "D"), label_size = 12, nrow=1)
b <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_stools_mock.png")
e <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_tissues_mock.png")
middle_panel <- plot_grid(b, e, labels = c("B", "E"), label_size = 12, ncol=2)
c <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_mock_stools.png")
f <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_mock_tissues.png")
bottom_panel <- plot_grid(c, f, labels = c("C", "F"), label_size = 12, nrow=1)

plot_grid(top_panel, middle_panel, bottom_panel, labels = NULL, label_size = 12, ncol=1, rel_heights = c(1.4, 1, 1))+
  ggsave("results/figures/figure_S3.pdf", width=6.875, height=7)+
  ggsave("submission/figure_S3.pdf", width=6.875, height=7)