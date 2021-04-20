source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_stool_PCoA_mock.png")
a_legend <- ggdraw() + draw_image("results/figures/5_days_PEG_pcoa_mock_legend.png")
a <- plot_grid(a, a_legend, ncol = 1, rel_heights = c(1, .2))
b <- ggdraw() + draw_image("results/figures/5_days_PEG_tissue_PCoA_mock.png")
b_legend <- ggdraw() + draw_image("results/figures/5_days_PEG_pcoa_mock_alpha_legend_tissues.png")
b <- plot_grid(b, b_legend, ncol = 1, rel_heights = c(1, .1))
top_panel <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, nrow=1)
c <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_stools_mock.png")
d <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_tissues_mock.png")
middle_panel <- plot_grid(c, d, labels = c("C", "D"), label_size = 12, ncol=2)
e <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_mock_stools.png")
f <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_mock_tissues.png")
bottom_panel <- plot_grid(e, f, labels = c("E", "F"), label_size = 12, nrow=1)
                       
plot_grid(top_panel, middle_panel, bottom_panel, labels = NULL, label_size = 12, ncol=1, rel_heights = c(1.4, 1, 1))+
  ggsave("results/figures/figure_S2.pdf", width=6.875, height=7)+
  ggsave("submission/figure_S2.pdf", width=6.875, height=7)
