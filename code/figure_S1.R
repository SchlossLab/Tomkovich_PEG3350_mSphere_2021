source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_cfu_WMR.png")
b <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_stools_WMR.png")

plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S1.pdf", width=6.875, height=7)+
  ggsave("submission/figure_S1.pdf", width=6.875, height=7)
