source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files
source("code/subset_analysis.R") #Read in PCoA subsets

#Define color scheme to match post_CDI_PEG Plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "CWM", "FRM", "RM")
color_labels <- c( "Clind.", "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350")


# Pull in diversity for alpha diversity analysis using post CDI PEG subset from defined in utilties.R
diversity_data_subset <- semi_join(diversity_data, post_cdi_PEG_metadata, by = c("unique_label")) #Only the samples that correspond to the post CDI PEG subset

#Statistical Analysis----
set.seed(19881117) #Match seed used in mothur analysis scripts

#Alpha Diversity Shannon Analysis----
shannon_post_cdi_peg <- diversity_data_subset %>%
  group_by(group, day) %>%
  mutate(median_shannon = median(shannon)) %>%
  ggplot(aes(x=group, y=shannon, colour=group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  #  scale_shape_manual(name=NULL,
  #                     values=shape_scheme,
  #                     breaks=shape_experiment,
  #                     labels=shape_experiment) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median_shannon, ymin = median_shannon), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(size=2, alpha=0.6, show.legend = TRUE) +
  labs(title=NULL,
       x=NULL,
       y="Shannon Diversity Index")+
  ylim(0, 3.5)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = c(.75,.25),
        text = element_text(size = 14), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank())
save_plot("results/figures/shannon_post_CDI_PEG.png", shannon_post_cdi_peg)

#Plot Shannon over time for all 30 days for post CDI PEG subset
shannon_post_cdi_peg_overtime_full <- plot_shannon_overtime(diversity_data_subset) +
  scale_x_continuous(breaks = c(-1:10, 15, 20, 25, 30),
                     limits = c(-2,35),
                     minor_breaks = c(-2.5:7.5)) +
  theme(legend.position = c(1, .15))
save_plot("results/figures/shannon_post_cdi_peg_overtime.png", shannon_post_cdi_peg_overtime_full)

#Plot Shannon over time for first 10 days for post CDI PEG subset
diversity_data_subset_10d <- diversity_data_subset %>% filter(day %in% c(-1:10))
shannon_post_cdi_peg_overtime_10d <- plot_shannon_overtime(diversity_data_subset_10d) +
  scale_x_continuous(breaks = c(-1:10),
                     limits = c(-2,17),
                     minor_breaks = c(-2.5:7.5)) +
  theme(legend.position = c(.92, .11))
#Haven;t finalized plot as of 10/30
#save_plot("results/figures/shannon_post_cdi_peg_overtime.png", shannon_post_cdi_peg_overtime_full)