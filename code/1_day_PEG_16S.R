source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match 1 Day Peg Plots
color_scheme <- c("#225ea8", "#238b45", "#88419d") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c( "1RM1", "C", "M1")
color_labels <- c( "1-day PEG 3350 + 1-day recovery", "Clind.", "1-day PEG 3350")

#Read in and create 1_Day_PEG metadata
metadata <- metadata %>%
  filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon")) %>% #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected
  filter(group == "C" & exp_num %in% c("M6")| #Only use C mice from this experiments. Allocated groups to figures based on paper outline.
           group == "1RM1" & exp_num %in% c("M6R")| #Had to differentiate experiment 6 from 6R in the metadata to create unique_mouse_id that wouldn't overlap for the M1 & 1RM1 mice that are both labeled with mouse_ids that are #s1-6
           group == "M1" & exp_num %in% c("M6"))%>%
  mutate(group=factor(group, levels=c("C", "1RM1", "M1")))  # Make sure group is treated as a factor

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


#Plot Shannon Diversity overtime
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
  

#Plot sobs for 1_day_PEG subset
  sobs_1_Day_PEG <- diversity_data_subset %>%
    group_by(group, day) %>%
    mutate(median_sobs = median(sobs)) %>%
    ggplot(x = day, y = sobs, colour =  group) +
    geom_point(mapping = aes(x = day, y = sobs, color = group, fill = group), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mapping = aes(x = day, y = median_sobs, color = group), alpha = 0.6, size = 1) +
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
  save_plot("results/figures/richness_1_Day_PEG.png", sobs_1_Day_PEG)
  
  
  
  
  