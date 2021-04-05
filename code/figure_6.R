source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Figure with Schematic, CFU over time, and shannon over time
a <- ggdraw() + draw_image("results/figures/ml_performance_genus.png")
b <- ggdraw() + draw_image("results/figures/ml_top_features_genus.png")
top_panel <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, nrow=1, rel_widths = c(1,2))
c <- ggdraw() + draw_image("results/figures/ml_d5_top10_genus.png")#Note scale could be removed

plot_grid(top_panel, c, labels = c("", "C"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_6.pdf", width=8, height=8)+
  ggsave("submission/figure_6.pdf", width=8, height=8)
