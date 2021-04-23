source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/ml_top_features_genus.png")
b <- ggdraw() + draw_image("results/figures/ml_top10_d5_genus.png")
c <- ggdraw() + draw_image("results/figures/ml_abund_5_genus_lineplot.png")

plot_grid(a, b, c, labels = c("A", "B", "C"), label_size = 12, ncol=1, rel_heights = c(2, 2, 1))+
  ggsave("results/figures/figure_6.pdf", width=6.8768, height=9)+
  ggsave("submission/figure_6.pdf", width=6.876, height=9)

