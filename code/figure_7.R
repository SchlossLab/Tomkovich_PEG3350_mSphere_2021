source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/summary_schematic.png")

plot_grid(a, labels = NULL, label_size = 12, ncol=1)+
  ggsave("results/figures/figure_7.pdf", width=6.875, height=4.5)+
  ggsave("submission/figure_7.pdf", width=6.875, height=4.5)

