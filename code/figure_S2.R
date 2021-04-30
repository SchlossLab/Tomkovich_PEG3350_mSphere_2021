source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_cecum_d6_10_genera.png")
b <- ggdraw() + draw_image("results/figures/5_days_PEG_genus_lineplot_peptostreptococcacea.png")
b <- plot_grid(b, NULL, labels = NULL)

plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol=1, rel_heights = c(1.4, 1, 1))+
  ggsave("results/figures/figure_S2.pdf", width=6.875, height=7)+
  ggsave("submission/figure_S2.pdf", width=6.875, height=7)

