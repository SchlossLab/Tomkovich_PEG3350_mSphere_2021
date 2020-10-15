source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match 1 Day Peg Plots
color_scheme <- c("#225ea8", "#238b45", "#88419d") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c( "1RM1", "C", "M1")
color_labels <- c( "1-day PEG 3350 + 1-day recovery", "Clind.", "1-day PEG 3350")

#Read in 1_Day_PEG metadata
metadata <- one_day_PEG_metadata


# Pull in diversity for alpha diversity analysis
diversity_data <- semi_join(diversity_data, metadata, by = c("unique_label")) #Only the samples that correspond to the 1_Day_PEG Subset

#Alpha Diversity Analysis - Shannon Plot
shannon_1_Day_PEG <- diversity_data %>%
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
save_plot("results/figures/shannon_1_Day_PEG.png", shannon_1_Day_PEG)

diversity_data_subset <- diversity_data %>%
  filter(day %in% c(-2, -1, 0, 1, 2, 4, 7))

Shannon_1_Day_PEG_Overtime <- diversity_data_subset %>%
  group_by(group, day) %>%
  mutate(median_shannon = median(shannon)) %>%
  ggplot(x = day, y = shannon, colour = group)+
  geom_point(mapping = aes(x = day, y = shannon, color = group, fill = group), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
  geom_line(mapping = aes(x = day, y = median_shannon, color = group), alpha = 0.6, size = 1) +
  scale_colour_manual(name=NULL,
                    values=color_scheme,
                    breaks=color_groups,
                    labels=color_labels) +
  scale_x_continuous(breaks = c(-2:7),
                     limits = c(-3, 8),
                     minor_breaks = c(-2.5:7.5))+
  theme_classic()+
  theme(legend.position = c(1,.25),
        text = element_text(size = 14), # Change font size for entire plot
        axis.ticks.x = element_blank())
  save_plot("results/figures/shannon_1_Day_PEG_Overtime.png", Shannon_1_Day_PEG_Overtime)
  
