source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_cfu_WMR_indiv.png")
b <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_stools_WMR_indiv.png")

plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol=1, rel_heights = c(1, 2))+
  ggsave("results/figures/figure_S1.pdf", width=6.875, height=9)+
  ggsave("submission/figure_S1.pdf", width=6.875, height=9)
