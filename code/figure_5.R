source("code/utilities.R") #Loads libraries, reads in metadata, functions
#One figure with schematic, CFU over time, and shannon over time (for stools only)? Check dimensions of Vendor paper for long wide version. plot grid for three panels that match up, and do another one? Vendor paper fig 1 can choose which panel goes were. 
#One figure with important 16S highlights: still pcoa, heatmap? Or choose a select couple days and do top genera/OTUs? 
#heat map main diffs between groups, keep it at genera level that don't recover in all PEG groups that remain colonized or have recovered by day 15
#Should I do a whole another fig for the tissue differences? Or should I combine with the 16S fig above?
#OTU12 over time probably in supplements bc it isnt as reliable
#Vendor paper Fig 6 makes fig legend a lot smaller
#Correlate c diff with microbiome changes
#Tropini 2018 says taxa that got effected by laxative treatments (bacteriodes, S24-7)

PEG_treatment_legend <- ggplot() +
  theme_void() +
  scale_x_continuous(limits = c(-10, 10))+
  scale_y_continuous(limits = c(-1.5, 1.5))+
  annotate("rect", xmin = -10, xmax = -9, ymin = -.75, ymax = .75, fill = "#88419d", alpha = .2)+ #shade to indicate PEG treatment in Clind + 1-day PEG group
  annotate("rect", xmin = -1.5, xmax = -.5, ymin = 0, ymax = .75, fill = "#f768a1", alpha = .2)+ #shade to indicate PEG treatment in Clind + 3-day recovery + 1-day PEG + FMT/PBS
  annotate("rect", xmin = -1.5, xmax = -.5, ymin = -.75, ymax = 0, fill = "#225ea8", alpha = .2)+
  annotate("text", x = -5.25, y = 0, label = "PEG treatment in Clind + 1-day PEG group", size = 3.78)+
  annotate("text", x = 5.3, y = 0, label = "PEG treatment in Clind + 3-day recovery + 1-day PEG (+ FMT/PBS)", size = 3.78)+
  xlab(NULL) 
save_plot("results/figures/post_CDI_PEG_treatment_legend.png", PEG_treatment_legend, base_height = 1, base_width = 9)
 
#Figure with Schematic, CFU over time, and shannon over time
a <- ggdraw() + draw_image("results/pictures/post_CDI_PEG_exp_schematic.png")
b <- ggdraw() + draw_image("results/figures/post_CDI_PEG_cfu.png")
c <- ggdraw() + draw_image("results/figures/post_CDI_PEG_shannon_stool.png")
d <- ggdraw() + draw_image("results/figures/post_CDI_PEG_richness_overtime_stool.png")
e <- ggdraw() + draw_image("results/figures/post_CDI_PEG_treatment_legend.png")
fig <- image_graph(width = 400, height = 400, res = 96)
plot_grid(a, e, b, c, d, labels = c("A", "", "B", "C", "D"), label_size = 12, ncol=1, rel_heights = c(.8,.2, 1,1,1))+
  ggsave("results/figures/figure_5.pdf", width=5.5, height=9)+
  ggsave("submission/figure_5.tiff", width=5.5, height=9, dpi = 600, device = "tiff", compression = "lzw", units = "in")
  
