source("code/utilities.R")

a <- ggdraw() + draw_image("results/pictures/1_Day_PEG_Schematic.png")
b <- ggdraw() + draw_image("results/figures/1_day_PEG_cfu.png")
c <- ggdraw() + draw_image("results/figures/1_day_PEG_otu_peptostreptococcaceae.png")
d <- ggdraw() + draw_image("results/figures/1_Day_PEG_shannon.png")
left_panel <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol = 1)
right_panel <- plot_grid(c, d, labels = c("C", "D"), label_size = 12, ncol = 1)
plot_grid(left_panel, right_panel) +
  ggsave("results/figures/figure_3.pdf", width=6.875, height=4.73)
