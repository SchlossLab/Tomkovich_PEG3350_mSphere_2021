source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/ml_abund_otu_muribaculum_lineplot.png")
b <- ggdraw() + draw_image("results/figures/ml_abund_otu_porphyromonadaceae_lineplot.png")
c <- ggdraw() + draw_image("results/figures/ml_abund_otu_lachnospiraceae_lineplot.png")

plot_grid(a, b, c, labels = c("A", "B", "C"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S3.pdf", width=6.875, height=6.5)+
  ggsave("submission/figure_S3.pdf", width=6.875, height=6.5)
