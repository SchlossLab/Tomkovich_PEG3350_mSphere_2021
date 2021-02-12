source("code/utilities.R") #Loads libraries, reads in metadata, functions
#One figure with schematic, CFU over time, and shannon over time (for stools only)? Check dimensions of Vendor paper for long wide version. plot grid for three panels that match up, and do another one? Vendor paper fig 1 can choose which panel goes were. 
#One figure with important 16S highlights: still pcoa, heatmap? Or choose a select couple days and do top genera/OTUs? 
#heat map main diffs between groups, keep it at genera level that don't recover in all PEG groups that remain colonized or have recovered by day 15
#Should I do a whole another fig for the tissue differences? Or should I combine with the 16S fig above?
#OTU12 over time probably in supplements bc it isnt as reliable
#Vendor paper Fig 6 makes fig legend a lot smaller
#Correlate c diff with microbiome changes
#Tropini 2018 says taxa that got effected by laxative treatments (bacteriodes, S24-7)

#Figure with Schematic, CFU over time, and shannon over time
a <- ggdraw() + draw_image("results/pictures/post_CDI_PEG_schematic.png")
b <- ggdraw() + draw_image("results/figures/post_CDI_PEG_cfu.png")
c <- ggdraw() + draw_image("results/figures/post_CDI_PEG_shannon_stool.png")
d <- ggdraw() + draw_image("results/figures/post_CDI_PEG_richness_overtime_stool.png")


plot_grid(a, b, c, d, labels = c("A", "B", "C", "D"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_4.pdf", width=5, height=7.5)


#Figure showing pcoa and genera heat map over time
d <- ggdraw() + draw_image("results/figures/post_CDI_PEG_stool_pcoa.png")
e <- ggdraw() + draw_image("results/figures/post_CDI_PEG_pcoa_legend.png")
pcoa <- plot_grid(d, e, label_size = 12, ncol = 1, rel_heights = c(3, 1))
heatmap <- ggdraw() + draw_image("results/figures/post_CDI_PEG_genus_heatmap_stools.png")

plot_grid(pcoa, heatmap, labels = c("A", "B"), ncol = 2, rel_heights = c(1, 3), rel_widths = c(.75, 1))+
  ggsave("results/figures/figure_4_16S.pdf", width = 10, height = 6.67)
