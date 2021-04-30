source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Figure with Schematic, CFU over time, and shannon over time
a <- ggdraw() + draw_image("results/figures/5_days_PEG_exp_schematic.png")
b <- ggdraw() + draw_image("results/figures/5_days_PEG_weight_median.png")
c <- ggdraw() + draw_image("results/figures/5_days_PEG_cfu.png")#Note scale could be removed

plot_grid(a, b, c, labels = c("A", "B", "C"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_1.pdf", width=6.875, height=8.375)+
  ggsave("submission/figure_1.pdf", width=6.875, height=8.375)
  
