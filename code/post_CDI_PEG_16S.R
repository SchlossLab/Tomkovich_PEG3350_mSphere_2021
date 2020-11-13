source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match post_CDI_PEG Plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8", "7f5f1e") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "CWM", "FRM", "RM", "FMT")
color_labels <- c( "Clind.", "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350", "FMT")

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Alpha Diversity Shannon Analysis----

# Pull in diversity for alpha diversity analysis using post CDI PEG subset from defined in utilties.R
diversity_data_subset <- post_cdi_PEG_subset(diversity_data) %>%
  add_row(diversity_data %>% filter(str_detect(unique_label, "FMT")) %>% #Also add the FMT gavage samples to this subset
            mutate(group = as.factor("FMT")))

shannon_post_cdi_peg <- diversity_data_subset %>%
  filter(group != "FMT") %>% #drop FMT from shannon
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
  ylim(0, 4)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = c(.5, .5), #Save Legend in center for later legend extraction
        text = element_text(size = 14), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank())
shannon_post_cdi_peg_no_legend <- shannon_post_cdi_peg +
  theme(legend.position = "none") #Remove legend
save_plot("results/figures/post_CDI_PEG_shannon.png", shannon_post_cdi_peg_no_legend) #Save Shannon plot without legend

#Extract Shannon Legend for all plots
legend_shannon_post_cdi_peg <- get_legend(shannon_post_cdi_peg) %>% as_ggplot()
save_plot("results/figures/post_CDI_PEG_shannon_legend.png", legend_shannon_post_cdi_peg)

#Plot Shannon over time for all 30 days for post CDI PEG subset
shannon_post_cdi_peg_overtime_full <- diversity_data_subset %>%
  filter(group != "FMT") %>% #drop FMTs
  plot_shannon_overtime() +
  scale_x_continuous(breaks = c(-1:10, 15, 20, 25, 30),
                     limits = c(-2,35),
                     minor_breaks = c(-2.5:7.5)) +
  scale_y_continuous(limits = c(0,4))+
  labs(x = "Day",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none") #Removing legend to save separately
save_plot("results/figures/post_CDI_PEG_shannon_overtime.png", shannon_post_cdi_peg_overtime_full) #Save full Shannon over time plot without legend


#Plot Shannon over time for first 10 days for post CDI PEG subset
diversity_data_subset_10d <- diversity_data_subset %>%
  filter(group != "FMT", #Drop FMTs
    day %in% c(-1:10))
shannon_post_cdi_peg_overtime_10d <- diversity_data_subset_10d %>%
  plot_shannon_overtime() +
  scale_x_continuous(breaks = c(-1:10),
                     limits = c(-2,11),
                     minor_breaks = c(-2.5:7.5)) +
  scale_y_continuous(limits = c(0,4)) +
  labs(x = "Day",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none")
save_plot("results/figures/shannon_post_cdi_peg_overtime_10d.png", shannon_post_cdi_peg_overtime_10d) #Save 10 day Shannon plot without legend





#Plot Stool + Tissue PCoA data----
#Pull post_CDI_PEG subset of PCoA data
pcoa_post_cdi_peg <- read_tsv("data/process/post_CDI_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  mutate(day = as.integer(day)) %>% #Day variable (transformed to integer to get rid of decimals on PCoA animation
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

#Pull axes from loadings file
pcoa_axes_post_cdi_PEG <- read_tsv("data/process/post_CDI_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_post_cdi_PEG %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_post_cdi_PEG %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

#PCoA plot and save the plot
pcoa_subset_plot <- plot_pcoa(pcoa_post_cdi_peg)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))
save_plot(filename = paste0("results/figures/post_CDI_PEG_pcoa.png"), pcoa_subset_plot, base_height = 7, base_width = 14)

#PCoA plot over time as a still and a as an animation
pcoa_plot_time <- plot_pcoa(pcoa_post_cdi_peg)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  theme( legend.position = "none")+ #remove legend
  facet_wrap(~ day)


pcoa_animated <- pcoa_post_cdi_peg %>%
  filter(group != "FMT") %>%
  plot_pcoa()+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Anotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  labs(title = 'Day: {frame_time}') + #Adds time variable to title
  transition_time(day)+  #Day variable used to cycle through time on animation
  shadow_mark() #Shows previous timepoints

# Implement better frames per second for animation
pcoa_gif <- animate(pcoa_animated, duration = , fps = 10,
                    res = 150, width = 20, height = 20, unit = "cm")
# Save as gif file
anim_save(animation = pcoa_gif, filename = 'results/post_CDI_PEG_pcoa_over_time.gif')

#Plot Stool Only PCoA data----
pcoa_post_cdi_peg_stool <- read_tsv("data/process/post_CDI_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  mutate(day = as.integer(day)) %>% #Day variable (transformed to integer to get rid of decimals on PCoA animation
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

#Pull axes from loadings file
pcoa_axes_post_cdi_PEG_stool <- read_tsv("data/process/post_CDI_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1_stool <- pcoa_axes_post_cdi_PEG_stool %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2_stool <- pcoa_axes_post_cdi_PEG_stool %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot_stool <- plot_pcoa(pcoa_post_cdi_peg_stool)+
  labs(x = paste("PCoA 1 (", axis1_stool, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2_stool,"%)", sep = ""))
save_plot(filename = paste0("results/figures/post_CDI_PEG_stool_pcoa.png"), pcoa_subset_plot_stool, base_height = 7, base_width = 14)

#Animate stool only PCoA over time
pcoa_animated_stool <- pcoa_post_cdi_peg_stool %>%
  filter(group != "FMT") %>%
  plot_pcoa()+
  labs(x = paste("PCoA 1 (", axis1_stool, "%)", sep = ""), #Anotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2_stool,"%)", sep = ""))+
  labs(title = 'Day: {frame_time}') + #Adds time variable to title
  transition_time(day)+  #Day variable used to cycle through time on animation
  shadow_mark() #Shows previous timepoints

# Implement better frames per second for animation
pcoa_gif_stool <- animate(pcoa_animated, duration = , fps = 10,
                    res = 150, width = 20, height = 20, unit = "cm")
# Save as gif file
anim_save(animation = pcoa_gif_stool, filename = 'results/post_CDI_PEG_stool_pcoa_over_time.gif')

#Plot Tissue Only PCoA data----
pcoa_post_cdi_peg_tissues <- read_tsv("data/process/post_CDI_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  mutate(day = as.integer(day)) %>% #Day variable (transformed to integer to get rid of decimals on PCoA animation
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff
#Pull axes from loadings file
pcoa_axes_post_cdi_PEG_tissues <- read_tsv("data/process/post_CDI_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1_tissues <- pcoa_axes_post_cdi_PEG_tissues %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2_tissues <- pcoa_axes_post_cdi_PEG_tissues %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot_tissues <- plot_pcoa(pcoa_post_cdi_peg_tissues)+
  labs(x = paste("PCoA 1 (", axis1_tissues, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2_tissues,"%)", sep = ""))
save_plot(filename = paste0("results/figures/post_CDI_PEG_tissue_pcoa.png"), pcoa_subset_plot_tissues, base_height = 7, base_width = 14)
