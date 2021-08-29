source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_dist_cecum_d6_10_genera.png")
b <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_peptostreptococcacea.png")
b <- plot_grid(b, NULL, labels = NULL, rel_widths = c(1.5, 1))

plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol=1, rel_heights = c(1.4, 1, 1))+
  ggsave("results/figures/figure_S2.pdf", width=6.875, height=7)+
  ggsave("submission/figure_S2.tiff", width=6.875, height=7, dpi = 600, device = "tiff", compression = "lzw", units = "in")

