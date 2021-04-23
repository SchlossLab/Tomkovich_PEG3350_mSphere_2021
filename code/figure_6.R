source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Vert version of Fig 4 16S
a <- ggdraw() + draw_image("results/figures/post_CDI_PEG_stool_pcoa.png", scale = 1.2)
b <- ggdraw() + draw_image("results/figures/post_CDI_PEG_pcoa_legend_vert.png")
c <- ggdraw() + draw_image("results/figures/post_CDI_PEG_genus_lineplot_fmt.png")
lineplot <- ggdraw() + draw_image("results/figures/post_CDI_PEG_genus_lineplot_stools.png")
pcoa1 <- plot_grid(a ,b, labels = c("A", ""), label_size = 12, nrow = 1, rel_heights = c(1,.5), rel_widths = c(1,.5))

plot_grid(pcoa1, c, lineplot, labels = c("", "B", "C"), label_size = 12, nrow = 3)+
  ggsave("results/figures/figure_6.pdf", width = 6.3, height = 9)+
  ggsave("submission/figure_6.pdf", width = 6.3, height = 9)
