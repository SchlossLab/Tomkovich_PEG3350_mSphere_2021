source("code/utilities.R") #Loads libraries, reads in metadata, functions

a <- ggdraw() + draw_image("results/figures/5_days_PEG_shannon_stools_mock.png")

plot_grid(a, labels = c("A"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S2.pdf", width=6.875, height=6)+
  ggsave("submission/figure_S2.pdf", width=6.875, height=6)
