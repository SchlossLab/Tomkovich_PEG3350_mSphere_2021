source("code/utilities.R")

a <- ggdraw() + draw_image("results/figures/1_Day_PEG_Schematic.png", scale = 1.45, x = .09)
b <- ggdraw() + draw_image("results/figures/1_day_PEG_cfu.png", scale = 1.00)
c <- ggdraw() + draw_image("results/figures/1_Day_PEG_PCoA.png", scale = 1)
legend <- ggdraw() + draw_image("results/figures/1_day_PEG_pcoa_legend.png", scale = 1)
c <- plot_grid(c, legend, nrow = 1, rel_widths = c(1, .1), scale = 1.18)
d <- ggdraw() + draw_image("results/figures/1_Day_PEG_shannon.png")
e <- ggdraw() + draw_image("results/figures/1_Day_PEG_genus_6_v2_baselinetoD1_lineplot.png")
top_panel <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, nrow = 1, rel_widths = c(.60, 1), rel_heights = c(1, 1))
middle_panel <- plot_grid(c, d, labels = c("C", "D"), label_size = 12, nrow = 1, rel_widths = c(.8, 1), rel_heights = c(1, 1))
top_panel <- plot_grid(top_panel, middle_panel, ncol = 1)
bottom_panel <- plot_grid(e, labels = c("E"), label_size = 12)
plot_grid(top_panel, bottom_panel, nrow = 2, rel_heights = c(1.1, .8)) +
  ggsave("results/figures/figure_4.pdf", width = 6.875, height = 6)+
  ggsave("submission/figure_4.tiff", width = 6.875, height = 6, ,  dpi = 600, device = "tiff", compression = "lzw", units = "in")

