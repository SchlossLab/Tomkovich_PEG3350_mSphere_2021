source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/ml_performance_genus.png")
b <- ggdraw() + draw_image("results/figures/ml_top_features_genus.png")
top_panel <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, nrow=1, rel_widths = c(1,1.5), rel_heights = c(2, 1))
c <- ggdraw() + draw_image("results/figures/ml_d5_top10_genus.png")#Note scale could be removed

plot_grid(top_panel, c, labels = c("", "C"), label_size = 12, ncol=1, rel_widths = c(1, 3), rel_heights = c(1, 1.5))+
  ggsave("results/figures/figure_6.pdf", width=8, height=7)+
  ggsave("submission/figure_6.pdf", width=8, height=7)

