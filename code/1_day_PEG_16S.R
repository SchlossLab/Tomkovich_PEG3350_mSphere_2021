source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files
source("code/subset_analysis.R") #Read in PCoA subsets


#Define color scheme to match 1 Day Peg Plots
color_scheme <- c("#225ea8", "#238b45", "#88419d") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c( "1RM1", "C", "M1")
color_labels <- c( "1-day PEG 3350 + 1-day recovery", "Clind.", "1-day PEG 3350")


# Pull in diversity for alpha diversity analysis using one day PEG subset from defined in utilties.R
diversity_data_subset <- semi_join(diversity_data, one_day_PEG_metadata, by = c("unique_label")) #Only the samples that correspond to the 1_Day_PEG Subset

#Statistical Analysis----
set.seed(19881117) #Match seed used in mothur analysis scripts

#Alpha Diversity Analysis - Shannon Plot 
shannon_1_Day_PEG <- one_day %>%
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
Shannon_1_Day_PEG_Overtime <- plot_shannon_overtime(diversity_data_subset) +
  scale_x_continuous(breaks = c(-2:7),
                     limits = c(-3,8),
                     minor_breaks = c(-2.5:7.5))
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
  
  
#Plot PCoA data----
  plot_pcoa <- function(df){
    ggplot(df, aes(x=axis1, y=axis2, color = group, alpha = day)) +
      geom_point(size=2) +
      scale_colour_manual(name=NULL,
                          values=color_scheme,
                          breaks=color_groups,
                          labels=color_labels)+
      #    scale_alpha_continuous(range = c(.3, 1),
      #                           breaks= c(2, 4, 6, 8, 10),
      #                           labels=c(2, 4, 6, 8, 10))+
      coord_fixed() +
      #    xlim(-0.4, 0.65)+
      #    ylim(-0.45, 0.6)+
      labs(x="PCoA 1",
           y="PCoA 2",
           alpha= "Day") +
      theme_classic()
  }
  
  #Pull 1_Day_PEG subset of PCoA data
  pcoa_1_day_PEG <- pcoa_data %>% 
    inner_join(one_day_PEG_metadata)
  
  #PCoA plot that combines the 2 experiments and save the plot----
  pcoa_plot <- plot_pcoa(pcoa_1_day_PEG)+
    labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
         y = paste("PCokA 2 (", all_axis2,"%)", sep = ""))
  save_plot(filename = paste0("results/figures/1_Day_PEG_PCoA.png"), pcoa_plot, base_height = 5, base_width = 4.5)
  
  pcoa_plot_time <- plot_pcoa(pcoa_1_day_PEG)+
    labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
         y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))+
    theme( legend.position = "none")+ #remove legend
    facet_wrap(~ day)
  
  