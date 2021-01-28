source("code/utilities.R") #Loads libraries, reads in metadata, functions
#One figure with schematic, CFU over time, and shannon over time (for stools only)? Check dimensions of Vendor paper for long wide version. plot grid for three panels that match up, and do another one? Vendor paper fig 1 can choose which panel goes were. 
#One figure with important 16S highlights: still pcoa, heatmap? diversity? Or choose a select couple days and do top genera/OTUs? 
#heat map main diffs between groups, keep it at genera level that don't recover in all PEG groups that remain colonized or have recovered by day 15
#Should I do a whole another fig for the tissue differences? Or should I combine with the 16S fig above?
#OTU12 over time probably in suppkements bc it isnt as reliable
#Vendor paper Fig 6 makes fig legend a lot smaller
#Correlate c diff with microbiome changes
#Tropini 2018 says taxa that got effected by laxative treatments (bacteriodes, S24-7)

a <- ggdraw() + draw_image("results/pictures/post_CDI_PEG_schematic.png")
b <- ggdraw() + draw_image("results/figures/post_CDI_PEG_cfu.png")
c <- ggdraw() + draw_image("results/figures/post_CDI_PEG_shannon_stool.png")


plot_grid(a, b, c, labels = c("A", "B", "C", ""), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_4.pdf", width=5, height=7.5)
