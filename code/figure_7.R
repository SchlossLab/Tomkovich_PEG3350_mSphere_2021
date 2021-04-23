source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/ml_top_features_genus.png")
b <- ggdraw() + draw_image("results/figures/ml_top10_d5_genus.png")
c <- ggdraw() + draw_image("results/figures/ml_abund_5_genus_lineplot.png")

plot_grid(a, b, c, labels = c("A", "B", "C"), label_size = 12, ncol=1, rel_heights = c(1, 2, 1))+
  ggsave("results/figures/figure_7.pdf", width=6.8, height=7.1)+
  ggsave("submission/figure_7.pdf", width=6.8, height=7.1)


