source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Figure with Schematic, CFU over time, and shannon over time
a <- ggdraw() + draw_image("results/figures/5_days_PEG_exp_schematic.png")
b <- ggdraw() + draw_image("results/figures/5_days_PEG_cfu.png")

plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_1.pdf", width=5, height=4.5)+
  ggsave("submission/figure_1.pdf", width=5, height=4.5)
  
