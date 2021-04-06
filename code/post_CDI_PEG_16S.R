source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match post_CDI_PEG Plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8", "7f5f1e") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "CWM", "FRM", "RM", "FMT")
color_labels <- c( "Clind.", "Clind. + 1-day PEG", "Clind. + 3-day recovery + 1-day PEG + FMT", "Clind. + 3-day recovery + 1-day PEG + PBS", "FMT")

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Alpha Diversity Analysis-------
# Pull in diversity for alpha diversity analysis using post CDI PEG subset from defined in utilties.R
#Diversity data for all days
diversity_data <- diversity_data %>%
  mutate(day = as.integer(day))
diversity_data_subset <- post_cdi_PEG_subset(diversity_data) %>%
  add_row(diversity_data %>% filter(str_detect(unique_label, "FMT")) %>% #Also add the FMT gavage samples to this subset
            mutate(group = as.factor("FMT")))

diversity_data_subset_10d <- diversity_data_subset %>%
  filter(group != "FMT", #Drop FMTs
         day %in% c(-1:10))

#Create subset of the post CDI PEG diversity data for just stool samples and tissues
diversity_stools <- subset_stool(diversity_data_subset)
diversity_tissues <- subset_tissue(diversity_data_subset)

#Experimental days to analyze with the Kruskal-Wallis test (timepoints with 16S data for at least 3 groups)
#Compare days with stool data for all groups
count_subset(diversity_stools) 
stool_test_days <- c(-1, 0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 15) #Days with stools in for all four groups

#Alpha Diversity Shannon Analysis----
#Function to perform Kruskal-Wallis test for differences in Shannon diversity index across groups on a particular day with Benjamini Hochberg correction
#Arguments: 
#diversity_subset <- subset (stools or tissue samples) of diversity_data to perform statistical test on
#timepoint = timepoints to assess differences between groups specific to the subset (stool or tissue)
#subset_name = label to append to results filename to indicate subset analyzed. Ex. stool, tissues, stool_mock, tissues_mock
kruskal_wallis_shannon <- function(diversity_subset, timepoints, subset_name){
  diversity_stats <- diversity_subset %>% 
    filter(day %in% timepoints) %>% 
    select(group, shannon, day) %>% #Embrace needed around diversity_measure to call a column in the dataframe
    group_by(day) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$shannon, g=as.factor(.x$group)) %>% tidy())) %>% 
    mutate(median = map(data, get_shannon_median_group)) %>% 
    unnest(c(model, median)) %>% 
    ungroup()  
  diversity_stats_adjust <- diversity_stats %>% 
    select(day, statistic, p.value, parameter, method, C, CWM, FRM, RM) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/post_CDI_PEG_shannon_stats_", subset_name, "_subset.tsv"))
}

#Test with shannon for stool subset
kw_shannon_stools <- kruskal_wallis_shannon(diversity_stools, stool_test_days, "stools")
sig_shannon_days_stools <- pull_sig_days(kw_shannon_stools)

#Plot Shannon
shannon_post_cdi_peg <- diversity_data_subset %>%
  filter(group != "FMT") %>% #drop FMT from shannon
  group_by(group, day) %>%
  mutate(median_shannon = median(shannon)) %>%
  ggplot(aes(x=group, y=shannon, colour=group, alpha = day))+
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

#Plot Shannon over time days -1 to 30 for post CDI PEG subset: all samples except FMTs
x_annotation <- sig_shannon_days_stools
y_position <- max(diversity_stools$shannon)+ 0.05
label <- kw_label(kw_shannon_stools)
shannon_post_cdi_peg_overtime_full <- diversity_data_subset %>%
  filter(group != "FMT") %>% #drop FMTs
  filter(day != "-15") %>% #Remove -15 timepoint
  plot_shannon_overtime() +
  scale_x_continuous(breaks = c(-1:10, 15, 20, 25, 30),
                     limits = c(-2,35), #removes day -15 here
                     minor_breaks = c(-2.5:10.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5)) +
  labs(x = "Day",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none") #Removing legend to save separately
save_plot("results/figures/post_CDI_PEG_shannon_overtime.png", shannon_post_cdi_peg_overtime_full) #Save full Shannon over time plot without legend


#Plot Shannon over time for first 10 days for post CDI PEG subset (all samples)
shannon_post_cdi_peg_overtime_10d <- diversity_data_subset_10d %>%
  plot_shannon_overtime() +
  scale_x_continuous(breaks = c(-1:10),
                     limits = c(-1.5,10.5),
                     minor_breaks = c(-1.5:10.5)) +
  labs(x = "Day",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none")
save_plot("results/figures/shannon_post_cdi_peg_overtime_10d.png", shannon_post_cdi_peg_overtime_10d) #Save 10 day Shannon plot without legend
#Warning because stat for day 15 annotation is removed when we limit x axis to day 1:10

#Plot Shannon over time for stools only
x_annotation <- sig_shannon_days_stools
y_position <- max(diversity_stools$shannon)+ 0.05
label <- kw_label(kw_shannon_stools)
shannon_post_cdi_peg_overtime_stool <- diversity_stools %>%
  filter(group != "FMT") %>% #drop FMTs
  plot_shannon_overtime() +
  scale_x_continuous(breaks = c(-1:10, 15, 20, 25, 30),
                     limits = c(-2,31), #removes day -15 here
                     minor_breaks = c(-1.5:10.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5)) +
  theme(legend.position = "none")+ #Removing legend to save separately
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = Inf, fill = "#88419d", alpha = .15)+ #shade to indicate PEG treatment in Clind + 1-day PEG group
  annotate("rect", xmin = 3, xmax = 4, ymin = 0, ymax = 2, fill = "#225ea8", alpha = .15)+ #shade to indicate PEG treatment in Clind + 3-day recovery + 1-day PEG + FMT/PBS
  annotate("rect", xmin = 3, xmax = 4, ymin = 2, ymax = Inf, fill = "#f768a1", alpha = .15)
 
save_plot("results/figures/post_CDI_PEG_shannon_stool.png", shannon_post_cdi_peg_overtime_stool, base_height = 4, base_width = 8.5, base_aspect_ratio = 2) #Save full Shannon over time plot without legend

#Alpha Diversity Richness (Sobs)----
#Function to perform Kruskal-Wallis test for differences in richness (sobs) across groups on a particular day with Benjamini Hochberg correction
#Arguments: 
#diversity_subset <- subset (stools or tissue samples) of diversity_data to perform statistical test on
#timepoint = timepoints to assess differences between groups specific to the subset (stool or tissue)
#subset_name = label to append to results filename to indicate subset analyzed. Ex. stool, tissues, stool_mock, tissues_mock
kruskal_wallis_richness <- function(diversity_subset, timepoints, subset_name){
  diversity_stats <- diversity_subset %>% 
    filter(day %in% timepoints) %>% 
    select(group, sobs, day) %>% #Embrace needed around diversity_measure to call a column in the dataframe
    group_by(day) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$sobs, g=as.factor(.x$group)) %>% tidy())) %>% 
    mutate(median = map(data, get_sobs_median_group)) %>% 
    unnest(c(model, median)) %>% 
    ungroup()  
  diversity_stats_adjust <- diversity_stats %>% 
    select(day, statistic, p.value, parameter, method, C, CWM, FRM, RM) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/post_CDI_PEG_richness_stats_", subset_name, "_subset.tsv"))
}
#Test with richness for stool subset
kw_richness_stools <- kruskal_wallis_richness(diversity_stools, stool_test_days, "stools")
sig_richness_days_stools <- pull_sig_days(kw_richness_stools)

#Plot Richness overtime 
sobs_post_CDI_PEG <- diversity_data_subset %>%
  filter(group != "FMT") %>% #Remove FMTs
  group_by(group, day) %>%
  mutate(median_sobs = median(sobs)) %>%
  ggplot(x = day, y = sobs, colour =  group) +
  geom_point(mapping = aes(x = day, y = sobs, color = group, fill = group), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
  geom_line(mapping = aes(x = day, y = median_sobs, color = group), alpha = 0.6, size = 1) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels) +
  scale_x_continuous(breaks = c(-15, -1:10, 15, 20, 25, 30),
                     limits = c(-15,35), 
                     minor_breaks = c(-15.5, -14.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5)) +
  theme_classic()+
  theme(legend.position = "none", #Remove legend
        panel.grid.minor.x = element_line(size = 0.4, color = "grey")) + #Add gray lines to clearly separate symbols by days
  labs(x = "Days Post-Infection",
       y = "Number of Observed OTUs")
save_plot("results/figures/post_CDI_PEG_richness_overtime.png", sobs_post_CDI_PEG)

#Richness oVertime 10 Day Version
sobs_post_CDI_PEG_10d <- diversity_data_subset_10d %>%
  filter(group != "FMT") %>% #Remove FMTs
  group_by(group, day) %>%
  mutate(median_sobs = median(sobs)) %>%
  ggplot(x = day, y = sobs, colour =  group) +
  geom_point(mapping = aes(x = day, y = sobs, color = group, fill = group), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
  geom_line(mapping = aes(x = day, y = median_sobs, color = group), alpha = 0.6, size = 1) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels) +
  scale_x_continuous(breaks = c(-1:10),
                     limits = c(-1.5,11),
                     minor_breaks = c(-1.5:10.5)) +
  theme_classic()+
  theme(legend.position = "none", #Remove legend
        panel.grid.minor.x = element_line(size = 0.4, color = "grey")) + #Add gray lines to clearly separate symbols by days
  labs(x = "Days Post-Infection",
       y = "Number of Observed OTUs")
save_plot("results/figures/post_CDI_PEG_richness_overtime_10d.png", sobs_post_CDI_PEG_10d)

#Richness over time for stools
x_annotation <- sig_richness_days_stools
y_position <- max(diversity_stools$sobs)+1.5
label <- kw_label(kw_richness_stools)
sobs_post_CDI_PEG_stool <- diversity_stools %>%
  filter(group != "FMT") %>% #Remove FMTs
  filter(day != "-15") %>% #Remove -15 timepoint
  plot_richness_overtime()+
  scale_x_continuous(breaks = c(-1:10, 15, 20, 25, 30),
                     limits = c(-2, 31),
                     minor_breaks = c(-1.5:10.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5)) +
 scale_y_continuous(limits = c(0,110))+
  theme_classic()+
  theme(legend.position = "none", #Remove legend
        panel.grid.minor.x = element_line(size = 0.4, color = "grey"))+ #Add gray lines to clearly separate symbols by days
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = Inf, fill = "#88419d", alpha = .15)+ #shade to indicate PEG treatment in Clind + 1-day PEG group
  annotate("rect", xmin = 3, xmax = 4, ymin = 0, ymax = 50, fill = "#225ea8", alpha = .15)+ #shade to indicate PEG treatment in Clind + 3-day recovery + 1-day PEG + FMT/PBS
  annotate("rect", xmin = 3, xmax = 4, ymin = 50, ymax = Inf, fill = "#f768a1", alpha = .15)
save_plot("results/figures/post_CDI_PEG_richness_overtime_stool.png", sobs_post_CDI_PEG_stool, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Plot Stool + Tissue PCoA data----
#Pull post_CDI_PEG subset of PCoA data
pcoa_post_cdi_peg <- read_tsv("data/process/post_CDI_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  mutate(day = as.integer(day)) %>% #Day variable (transformed to integer to get rid of decimals on PCoA animation
  filter(!is.na(axis1)) %>% #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff
  filter(day > -3| is.na(day))#limit to experimental time frame & FMT samples

#Pull axes from loadings file
pcoa_axes_post_cdi_PEG <- read_tsv("data/process/post_CDI_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_post_cdi_PEG %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_post_cdi_PEG %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

#PCoA plot and save the plot
pcoa_subset_plot <- plot_pcoa(pcoa_post_cdi_peg)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))
save_plot(filename = paste0("results/figures/post_CDI_PEG_pcoa.png"), pcoa_subset_plot, base_height = 5, base_width = 5)

#Create stand alone legend
group_legend <- pcoa_post_cdi_peg  %>%
  ggplot(aes(x = axis1, y = axis2, color = group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  geom_point()+ theme_classic()+
  guides(color = guide_legend(ncol = 2))
group_legend <- get_legend(group_legend)
save_plot("results/figures/post_CDI_PEG_pcoa_legend.png", group_legend, base_height = .8, base_width = 6)
#Create Standalone legend, vertical
group_legend_vert <- pcoa_post_cdi_peg  %>%
  ggplot(aes(x = axis1, y = axis2, color = group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels= str_wrap(color_labels, 25))+
  geom_point()+ theme_classic()
group_legend_vert <- get_legend(group_legend_vert)
save_plot("results/figures/post_CDI_PEG_pcoa_legend_vert.png", group_legend_vert, base_height = 1.5, base_width = 2)

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
  filter(!is.na(axis1)) %>% #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff
  filter(day > -3 | is.na(day))#limit to experimental time frame & FMT samples

#Pull axes from loadings file
pcoa_axes_post_cdi_PEG_stool <- read_tsv("data/process/post_CDI_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1_stool <- pcoa_axes_post_cdi_PEG_stool %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2_stool <- pcoa_axes_post_cdi_PEG_stool %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

#Import current R-squared values and P values for top contributing to PCoA
#permanova_results = read_tsv("data/process/post_CDI_PEG_permanova_stools.tsv")
#r_sq_ordered_perm = permanova_results %>% arrange(desc(R2))
# Select top 2 contributors (3 to include the total header at the top where R2=1)
#top_2_contrib = top_n(r_sq_ordered_perm, 3, R2) #top two R2 for this subset is group:unique_cage_no instead of group, unlike 1 day subset
#top_R2 = signif(top_2_contrib[2,6], 4)
#top_R2_title = top_2_contrib[2,1]
#sec_R2 = signif(top_2_contrib[3,6], 4)
#sec_R2_title = "Group: Cage" #Replace top_2_contrib[3,1] (originally written as group:unique_cage_no) with "Group: Cage" 
#make titles their own variable nad have them get edited into upper case separately?

pcoa_subset_plot_stool <- plot_pcoa(pcoa_post_cdi_peg_stool)+
  labs(x = paste("PCoA 1 (", axis1_stool, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2_stool,"%)", sep = ""))
  #annotate("text", x = -.6, y = .41, label = paste(str_to_title(top_R2_title)), size = 3.2) +
  #annotate("text", x = -.45, y = .41, label = paste(str_to_title(sec_R2_title)), size = 3.2) +
  #annotate("text", x = -.6, y = .375, label = paste("italic(R)^2: ", top_R2, sep = "" ), parse = TRUE, size = 3.2) +
  #annotate("text", x = -.45, y = .375, label =  paste("italic(R)^2: ", sec_R2, sep = "" ), parse = TRUE, size = 3.2) +
  #annotate("text", x = -.6, y = .33, label = "italic(P) < 0.05", parse = TRUE, size = 3.2) +
  #annotate("text", x = -.45, y = .33, label = "italic(P) < 0.05", parse = TRUE, size = 3.2)
save_plot(filename = paste0("results/figures/post_CDI_PEG_stool_pcoa.png"), pcoa_subset_plot_stool, base_height = 6, base_width = 6)

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
pcoa_gif_stool <- animate(pcoa_animated_stool, duration = , fps = 10,
                    res = 150, width = 20, height = 20, unit = "cm")
# Save as gif file
anim_save(animation = pcoa_gif_stool, filename = 'results/post_CDI_PEG_stool_pcoa_over_time.gif')

#Plot Tissue Only PCoA data----
pcoa_post_cdi_peg_tissues <- read_tsv("data/process/post_CDI_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  mutate(day = as.integer(day)) %>% #Day variable (transformed to integer to get rid of decimals on PCoA animation
  filter(!is.na(axis1)) %>% #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff
  filter(day > -3)#limit to experimental time frame
#Pull axes from loadings file
pcoa_axes_post_cdi_PEG_tissues <- read_tsv("data/process/post_CDI_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1_tissues <- pcoa_axes_post_cdi_PEG_tissues %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2_tissues <- pcoa_axes_post_cdi_PEG_tissues %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot_tissues <- pcoa_post_cdi_peg_tissues %>% 
  mutate(`Mouse Number` = as.factor(mouse_id)) %>%
  ggplot(aes(x=axis1, y=axis2, color = `Mouse Number`, shape = sample_type)) + #Did not specitfy day since all tissues are from Day 30
  geom_point(size=2) +
  coord_fixed() + 
  theme_classic()+
  labs(x = paste("PCoA 1 (", axis1_tissues, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2_tissues,"%)", sep = ""),
       alpha = "Day",
       shape = "Tissue Type")
save_plot(filename = paste0("results/figures/post_CDI_PEG_tissue_pcoa.png"), pcoa_subset_plot_tissues, base_height = 5, base_width = 5)

#OTU Analysis------
#Function to test at the otu level:
agg_otu_data_subset <- post_cdi_PEG_subset(agg_otu_data) %>% 
  filter(sample_type =="stool") #Exclude the other sample types and just perform test on the stools

agg_otu_data_tissues <- post_cdi_PEG_subset(agg_otu_data) %>% 
  filter(!sample_type =="stool") %>% 
  mutate(sample_type = fct_relevel(sample_type, "cecum", "proximal_colon", "distal_colon")) #Specify order of sample types

kruskal_wallis_otu <- function(timepoint){
  otu_stats <- agg_otu_data_subset %>%
    filter(day == timepoint) %>%
    select(group, otu, agg_rel_abund) %>%
    group_by(otu) %>%
    nest() %>%
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$group)) %>% tidy())) %>%
    mutate(median = map(data, get_rel_abund_median_group)) %>%
    unnest(c(model, median)) %>%
    ungroup()
  #Adjust p-values for testing multiple OTUs
  otu_stats_adjust <- otu_stats %>%
    select(otu, statistic, p.value, parameter, method, "C", "CWM", "FRM", "RM") %>%
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
    arrange(p.value.adj) %>%
    write_tsv(path = paste0("data/process/post_CDI_PEG_otu_stats_day_", timepoint, ".tsv"))
}

#Create empty data frame to combine stat dataframes for all days that were tested
kw_otu_stools <- data.frame(otu=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                            C =double(),CWM =double(),FRM =double(),RM=double(),
                            p.value.adj=double(),day=double())


## Perform kruskal wallis tests at the otu level for all days of the experiment that were sequenced----
#Stool samples
for (d in c(-1, 0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 15)){
  kruskal_wallis_otu(d)
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/post_CDI_PEG_otu_stats_day_", d, ".tsv"))%>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_otu_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, otu))
  kw_otu_stools <- add_row(kw_otu_stools, stats)  #combine all the dataframes together
}

#Shared significant OTUs across from Day 3, 5, 6, 8, 10
shared_sig_otus_D3toD10 <- intersect_all(sig_otu_day3, sig_otu_day5, sig_otu_day6, sig_otu_day8, sig_otu_day10)
view(shared_sig_otus_D3toD10)
print(shared_sig_otus_D3toD10)

#[1] "Lachnospiraceae (OTU 4)"     "Lachnospiraceae (OTU 11)"    "Lachnospiraceae (OTU 33)"   
#[4] "Lachnospiraceae (OTU 24)"    "Lachnospiraceae (OTU 31)"    "Lachnospiraceae (OTU 30)"   
#[7] "Porphyromonadaceae (OTU 14)" "Porphyromonadaceae (OTU 44)"

#Transform day into a factor label to plot bacteria over time
agg_otu_data_subset <- agg_otu_data_subset %>% 
  mutate(day = fct_relevel(day, "-15", "-1", "0", "1", "2", "3", "4", "5", "6", "7", 
                           "8", "9", "10", "15", "20", "25", "30"))

#Examine C. difficile OTU over time----
peptostrep_stools <- otu_over_time("Peptostreptococcaceae (OTU 12)", agg_otu_data_subset)+
  scale_x_discrete(breaks = c(-1:10, 15, 20, 25, 30), labels = c(-1:10, 15, 20, 25, 30)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/post_CDI_PEG_otu_peptostreptococcaceae.png", peptostrep_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
peptostrep_stools_v2 <- otu_over_time("Peptostreptococcaceae (OTU 12)", agg_otu_data_subset)+
  scale_x_discrete(breaks = c(-15, -1:10, 15, 20, 25, 30), labels = c(-15, -1:10, 15, 20, 25, 30)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/post_CDI_PEG_otu_peptostreptococcaceae_dn15.png", peptostrep_stools_v2, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
#C. diff OTU across sample types on day 30
peptostrep_tissues <- otu_gi_distrib("Peptostreptococcaceae (OTU 12)", agg_otu_data_tissues, "30", "CWM")
save_plot(filename = "results/figures/post_CDI_PEG_otu_peptostreptococcaceae_CWM_tissues.png", peptostrep_tissues, base_height = 4, base_width = 6)

#Heatmap of significant OTUs ranked by-------
#Rank OTUs by adjusted p-value
hm_sig_otus_p_adj <- kw_otu_stools %>% 
  filter(p.value.adj < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(otu) %>% 
  slice_head(n = 25) %>% 
  pull(otu)

hm_stool_days <- diversity_stools %>% distinct(day) %>% 
  filter(!day == "-15") %>% pull(day)
facet_labels <- c("Clind.", "Clind. + 1-day PEG", "Clind. + 3-day recovery + 1-day PEG + FMT", "Clind. + 3-day recovery + 1-day PEG +PBS") #Create descriptive labels for facets
names(facet_labels) <- c("C", "CWM", "FRM", "RM") #values that correspond to group, which is the variable we're faceting by
hm_stool <- hm_plot_otus(agg_otu_data_subset, hm_sig_otus_p_adj, hm_stool_days)+
  scale_x_discrete(breaks = c(-1:10, 15, 20, 25, 30), labels = c(-1:10, 15, 20, 25, 30)) 
save_plot(filename = "results/figures/post_CDI_PEG_otus_heatmap_stools.png", hm_stool, base_height = 14, base_width = 15)

#Plot heat map of OTUs that were significant in the 5 day subset
hm_5_day_otus <- c("Phenylobacterium (OTU 332)", "Lachnospiraceae (OTU 33)", "Ruminococcaceae (OTU 37)",
                   "Ruminococcaceae (OTU 98)", "Ruminococcaceae (OTU 65)", "Clostridium (OTU 51)",
                   "Bacteroides (OTU 1)", "Lactobacillus (OTU 13)", "Blautia (OTU 19)", 
                   "Ruminococcaceae (OTU 50)", "Ruminococcaceae (OTU 54)", "Ruminococcaceae (OTU 92)",
                   "Bifidobacterium (OTU 28)", "Oscillibacter (OTU 45)", "Lachnospiraceae (OTU 30)", 
                   "Lachnospiraceae (OTU 29)", "Lachnospiraceae (OTU 31)", "Lactobacillus (OTU 23)", 
                   "Lachnospiraceae (OTU 4)", "Peptostreptococcaceae (OTU 12)", "Lachnospiraceae (OTU 32)", 
                   "Enterobacteriaceae (OTU 2)", "Clostridium (OTU 22)", "Lactobacillus (OTU 9)", "Lachnospiraceae (OTU 16)")
hm_stool_5_day_otus <- hm_plot_otus(agg_otu_data_subset, hm_5_day_otus, hm_stool_days)+
  scale_x_discrete(breaks = c(-1:10, 15, 20, 25, 30), labels = c(-1:10, 15, 20, 25, 30)) 
save_plot(filename = "results/figures/post_CDI_PEG_otus_heatmap_stools_5_day_otus.png", hm_stool_5_day_otus, base_height = 14, base_width = 15)

#Plot heatmaps of the tissue samples (only collected on day 30)
#Only collected tissues from CWM group: "Clind + 1-day PEG 3350"
hm_tissues_days <- 30
facet_labels <- c("Cecum", "Proximal colon", "Distal colon") #Create descriptive labels for facets
names(facet_labels) <- c("cecum", "proximal_colon", "distal_colon") #values that correspond to group, which is the variable we're faceting by
hm_tissues <- hm_plot_tissues(agg_otu_data_tissues, hm_sig_otus_p_adj, hm_tissues_days)+
  scale_x_discrete(breaks = c(30), labels = c(30)) 
save_plot(filename = "results/figures/post_CDI_PEG_otus_heatmap_tissues.png", hm_tissues, base_height = 10, base_width = 8)
hm_tissues_5_day_otus <- hm_plot_tissues(agg_otu_data_tissues, hm_5_day_otus, hm_tissues_days)+
  scale_x_discrete(breaks = c(30), labels = c(30))
save_plot(filename = "results/figures/post_CDI_PEG_otus_heatmap_tissues_5_day_otus.png", hm_tissues_5_day_otus, base_height = 10, base_width = 8)

#Genus Analysis----
agg_genus_data_subset <- post_cdi_PEG_subset(agg_genus_data) %>%
  filter(sample_type =="stool") #Exclude the other sample types and just perform test on the stools

agg_genus_data_tissues <- post_cdi_PEG_subset(agg_genus_data) %>%
  filter(!sample_type =="stool") %>% #Exclude all stool samples and just perform test on tissues
  mutate(sample_type = fct_relevel(sample_type, "cecum", "proximal_colon", "distal_colon")) #Specify order of sample types

kruskal_wallis_genus <- function(timepoint){
  genus_stats <- agg_genus_data_subset %>%
    filter(day == timepoint) %>%
    select(group, genus, agg_rel_abund) %>%
    group_by(genus) %>%
    nest() %>%
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$group)) %>% tidy())) %>%
    mutate(median = map(data, get_rel_abund_median_group)) %>%
    unnest(c(model, median)) %>%
    ungroup()
  #Adjust p-values for testing multiple OTUs
  genus_stats_adjust <- genus_stats %>%
    select(genus, statistic, p.value, parameter, method, "C", "CWM", "FRM", "RM") %>%
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
    arrange(p.value.adj) %>%
    write_tsv(path = paste0("data/process/post_CDI_PEG_genus_stats_day_", timepoint, ".tsv"))
}

#Create empty data frame to combine stat dataframes for all days that were tested
kw_genus_stools <- data.frame(genus=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                            C =double(),CWM =double(),FRM =double(),RM=double(),
                            p.value.adj=double(),day=double())


## Perform kruskal wallis tests at the genus level for all days of the experiment that were sequenced----
#Stool samples
for (d in c(-1, 0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 15)){
  kruskal_wallis_genus(d)
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/post_CDI_PEG_genus_stats_day_", d, ".tsv"))%>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_genus_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, genus))
  kw_genus_stools <- add_row(kw_genus_stools, stats)  #combine all the dataframes together
}

#Shared significant genera across from Day 2, 3, 5, 6, 7, 8, 10, 15
shared_sig_genus_D2toD15 <- intersect_all(sig_genus_day2, sig_genus_day3, sig_genus_day5, sig_genus_day6, sig_genus_day7, sig_genus_day8, sig_genus_day9, sig_genus_day10, sig_genus_day15)
view(shared_sig_genus_D2toD15)
print(shared_sig_genus_D2toD15)
#"Ruminococcaceae Unclassified" "Akkermansia" 

##Shared significant genera across from Day 3, 5, 6, 8, 10 (same days as for OTU)
shared_sig_genus_D3toD10 <- intersect_all(sig_genus_day3, sig_genus_day5, sig_genus_day6, sig_genus_day8, sig_genus_day9, sig_genus_day10)
View(shared_sig_genus_D3toD10)
print(shared_sig_genus_D3toD10)
#[1] "Ruminococcaceae Unclassified"    "Lachnospiraceae Unclassified"    "Porphyromonadaceae Unclassified"
#[4] "Oscillibacter"                   "Akkermansia"                     "Clostridiales Unclassified"  

pull_significant_taxa(kw_genus_stools, genus) #134 sig taxa across all timepoints

# Perform pairwise Wilcoxan rank sum tests for genera that were significantly different across groups on a series of days----
pairwise_day_genus <- function(timepoint, sig_genus_dayX){
  genus_stats <- post_cdi_PEG_subset(agg_genus_data) %>% 
    filter(day == timepoint) %>%
    select(group, genus, agg_rel_abund) %>% 
    group_by(genus) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$group)) %>% tidy())) %>% 
    mutate(median = map(data, get_rel_abund_median_group)) %>% 
    unnest(c(model, median)) %>% 
    ungroup()
  pairwise_stats <- genus_stats %>% 
    filter(genus %in% sig_genus_dayX) %>% 
    group_by(genus) %>% 
    mutate(model=map(data, ~pairwise.wilcox.test(x=.x$agg_rel_abund, g=as.factor(.x$group), p.adjust.method="BH") %>% 
                       tidy() %>% 
                       mutate(compare=paste(group1, group2, sep="-")) %>% 
                       select(-group1, -group2) %>% 
                       pivot_wider(names_from=compare, values_from=p.value)
    )
    ) %>% 
    unnest(model) %>% 
    select(-data, -parameter, -statistic, -p.value) %>% #Get rid of p.value since it's the unadjusted version
    write_tsv(path = paste0("data/process/post_CDI_PEG_genus_stats_day_", timepoint, "_sig.tsv"))
  #Format pairwise stats to use with ggpubr package
  plot_format_stats <- pairwise_stats %>% 
    select(-method, -C, -CWM, -FRM, -RM) %>% 
    group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
    lapply(tidy_pairwise_genus) %>% 
    bind_rows() %>% 
    arrange(p.adj) %>% #Arrange by adjusted p value column 
    mutate(day = timepoint)
  return(plot_format_stats)  
}

#Day 2, 3, 5, 6, 7, 8, 10, 15
genus_day2_stats <- pairwise_day_genus(2, sig_genus_day2)
genus_day3_stats <- pairwise_day_genus(3, sig_genus_day3)
genus_day5_stats <- pairwise_day_genus(5, sig_genus_day5)
genus_day6_stats <- pairwise_day_genus(6, sig_genus_day6)
genus_day7_stats <- pairwise_day_genus(7, sig_genus_day7)
genus_day8_stats <- pairwise_day_genus(8, sig_genus_day8)
genus_day10_stats <- pairwise_day_genus(10, sig_genus_day10)
genus_day15_stats <- pairwise_day_genus(15, sig_genus_day15)
genus_pairwise_stools <- rbind(genus_day2_stats, genus_day3_stats, genus_day5_stats, genus_day6_stats, genus_day7_stats, genus_day8_stats, genus_day10_stats, genus_day15_stats)
genus_pairwise_stools_5plusdpi <- rbind(genus_day5_stats, genus_day6_stats, genus_day7_stats, genus_day8_stats, genus_day10_stats, genus_day15_stats)

#Heatmap of significant genera ranked by
#Rank genera by adjusted p-value
hm_sig_genera_p_adj <- kw_genus_stools %>% 
  filter(p.value.adj < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(genus) %>% 
  slice_head(n = 25) %>% 
  pull(genus)

hm_stool_days <- diversity_stools %>% distinct(day) %>% # Redundant if whole script has already been run, but if not, requires lines 18 to 27 to have already run 
  filter(!day == "-15") %>% pull(day)
facet_labels <- c("Clind.", "Clind. + 1-day PEG", "Clind. + 3-day recovery + 1-day PEG+ FMT", "Clind. + 3-day recovery + 1-day PEG + PBS") #Create descriptive labels for facets
names(facet_labels) <- c("C", "CWM", "FRM", "RM") #values that correspond to group, which is the variable we're faceting by
agg_genus_data_subset_hm <- agg_genus_data_subset %>% filter(!(day %in% c(4, 20, 25))) #drop day 20 and 25 (only data for one group)
hm_stool <- hm_plot_genus(agg_genus_data_subset_hm, hm_sig_genera_p_adj, hm_stool_days)+
  scale_x_discrete(breaks = c(-1:10, 15, 30), labels = c(-1:10, 15, 30)) 
save_plot(filename = "results/figures/post_CDI_PEG_genus_heatmap_stools.png", hm_stool, base_height = 7, base_width = 10)

#Plot heatmap significant genera over time facet by genus
facet_labels <- c("Peptostreptococcaceae Unclassified", "Clostridiales Unclassified", "Oscillibacter", "Bacteroides", "Acetatifactor", "Akkermansia") 
names(facet_labels) <- c("Peptostreptococcaceae Unclassified", "Clostridiales  Unclassified", "Oscillibacter", "Bacteroides", "Acetatifactor", "Akkermansia")
exp_groups <- c("RM", "FRM", "CWM", "C") #Arrange this way to match pcoa legend
exp_group_labels <- c("Clind. + 3-day recovery + 1-day PEG","Clind. + 3-day recovery + 1-day PEG + FMT", "Clind. + 1-day PEG + PBS", "Clind.")

hm_genera_facet <- hm_plot_genus_facet(agg_genus_data_subset_hm, exp_groups, facet_labels, hm_stool_days, exp_group_labels)+
  scale_x_discrete(breaks = c(-1:10, 15, 30), labels = c(-1:10, 15, 30))
save_plot(filename = "results/figures/post_CDI_PEG_genus_heatmap_facet.png", hm_genera_facet, base_height = 7, base_width = 15)

#Plot alt heatmap sig genera over time facet by genus
#pull sig genera from:
#Pairwise comparison between Clind. + 3-day recovery + 1-day PEG 3350 with and without FMT 

sig_genera_pw_5dpi <- genus_pairwise_stools_5plusdpi %>% filter(group1 %in% c( "FRM", "RM") & group2 %in% c("FRM", "RM")) %>% 
  filter(p.adj < .05) %>%
  arrange(p.adj)  %>%
  filter(duplicated(genus)) %>%#pull sig genera over mulitiple days
  pull(genus)
sig_genera_pw_5dpi #Probably do not make much difference with cdi clearance

#Sig genera between groups for a day
sig_genera_kw_multi_day <- kw_genus_stools %>% filter(day >= 5, p.value.adj  < .05) %>% 
  filter(duplicated(genus),
        !(genus %in% sig_genera_pw_5dpi)) %>% #drop genera that were sig diff between PBS/FMT
 count(genus)%>% arrange(desc(n)) %>% #put genera sig over most days at top
 filter(genus != "Unclassified") %>%
  top_n(6) %>%  #select top 6 for heatmap
 pull(genus)
#Create heat map
facet_labels_alt <- sig_genera_kw_multi_day #top 6 from exploratory below
hm_genera_facet_alt <- hm_plot_genus_facet(agg_genus_data_subset_hm, exp_groups, facet_labels_alt, hm_stool_days, exp_group_labels)+
  scale_x_discrete(breaks = c(-1:10, 15, 30), labels = c(-1:10, 15, 30))
save_plot(filename = "results/figures/post_CDI_PEG_genus_heatmap_facet_alt.png", hm_genera_facet_alt, base_height = 7, base_width = 15)

#Plot genera heatmaps of the tissue samples (only collected on day 30)
#Only collected tissues from CWM group: "Clind + 1-day PEG 3350"
hm_tissues_days <- 30
facet_labels <- c("Cecum", "Proximal colon", "Distal colon") #Create descriptive labels for facets
names(facet_labels) <- c("cecum", "proximal_colon", "distal_colon") #values that correspond to group, which is the variable we're faceting by
hm_tissues_genera <- hm_plot_tissues_genera(agg_genus_data_tissues, hm_sig_genera_p_adj, hm_tissues_days)+
  scale_x_discrete(breaks = c(30), labels = c(30)) 
save_plot(filename = "results/figures/post_CDI_PEG_genera_heatmap_tissues.png", hm_tissues_genera, base_height = 10, base_width = 8)

#Exploratory--------

genus_day3_stats %>% 
  filter(group1 %in% c( "FRM", "RM") & group2 %in% c("FRM", "RM")) %>%
  filter(p.adj < .05) %>% 
  pull(genus) #[1] "Clostridium XlVb"             "Lachnospiraceae Unclassified"

genus_day5_stats %>% 
  filter(group1 %in% c( "FRM", "RM") & group2 %in% c("FRM", "RM")) %>%
  filter(p.adj < .05) %>% 
  pull(genus) #[1] "Porphyromonadaceae Unclassified"

#Select Genera to show what's different across all groups
pairwise_genus_rank <- genus_day5_stats %>% filter(p.adj < .05) %>% 
  count(genus) %>% arrange(desc(n)) %>% #Rank genera according to # of groups with sig diff
  top_n(6) %>% #select top 6 for over time plots
  pull(genus)
kw_genus_stools %>% #Check to see if there were any genera that were sig in the kw test that were not sig in the pairwise comparisons
  filter(day == 5,p.value.adj < .05,
         !(genus %in% pairwise_genus_rank)) %>% pull(genus) #Only [1] "Turicibacter"
#Plot alt line plots faceted by genus
agg_genus_data_subset_hm <- agg_genus_data_subset %>% filter(!(day %in% c(4, 20, 25))) #drop day 20 and 25 (only data for one group)
exp_groups <- c("RM", "FRM", "CWM", "C") #Arrange this way to match pcoa legend
exp_group_labels <- c("Clind. + 3-day recovery + 1-day PEG + PBS","Clind. + 3-day recovery + 1-day PEG + FMT", "Clind. + 1-day PEG", "Clind.")
line_plot_stool_days <- diversity_stools %>% distinct(day) %>% 
  filter(!(day %in% c(-15, 30))) %>% pull(day)
genus_line_plot_facet <- agg_genus_data_subset_hm %>% 
  mutate(group = fct_relevel(group, exp_groups)) %>% #Specify the order of the groups
  filter(genus %in% pairwise_genus_rank) %>%
  filter(day %in% line_plot_stool_days) %>% 
  mutate(day = as.numeric(day)) %>% 
  group_by(group, genus, day) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  ggplot()+
  geom_line(aes(x = day, y=median, color=group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
#  scale_x_continuous(limits = c(-1,15), breaks = c(-1:10, 15), labels = c(-1:10, 15))+
#  scale_y_continuous(trans = "log10", limits = c(1/10900, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  #geom_hline(yintercept=1/1000, color="gray")+
  labs(title=NULL,
       x="Days Post-Infection",
       y="Relative Abundance (%)")+
  facet_wrap(~genus, nrow = 2, labeller = label_wrap_gen(width = 10), scales = "free")+
  theme_classic()+
  theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
        strip.text = element_text(face = "italic"),
        plot.title = element_markdown(hjust = 0.5), #Have only the genera names show up as italics
        text = element_text(size = 16),
        legend.position = "bottom")
save_plot("results/figures/post_CDI_PEG_genus_line_plot_facet.png", genus_line_plot_facet)

#List of top 10 significant genera in stool samples that are significant over the most timepoints
top_10_sig_genus <- kw_genus_stools %>% 
  filter(p.value.adj < 0.05) %>% #Select only significant p-values
  group_by(genus) %>% 
  tally() %>% 
  filter(n > 1) %>% #select genera that vary over multiple timepoints = n > 1
  arrange(desc(n)) %>% #Rank by how many time each genera shows up (
  filter(!genus == "Unclassified") %>% #Remove this genus since it's not informative and could contain multiple unclassified genera
  head(10) %>% 
  pull(genus)

#Line plot of top 6 significant genera in stool samples that uses function in code/utilities.R
facet_labels <- top_10_sig_genus[1:6] #Pick just the top 6
names(facet_labels) <- top_10_sig_genus[1:6] #Pick just the top 6
line_plot_stool_days <- diversity_stools %>% distinct(day) %>% 
  filter(!(day %in% c(-15, 30, 20, 25))) %>% pull(day)
lp_stool <- line_plot_genus(agg_genus_data_subset_hm, top_10_sig_genus[1:6], line_plot_stool_days, "solid")+
  scale_x_continuous(limits = c(-1,15), breaks = c(-1:10, 15), labels = c(-1:10, 15)) #Rewrite over -1:10 scale
save_plot(filename = "results/figures/post_CDI_PEG_genus_lineplot_stools.png", lp_stool, base_height = 5, base_width = 8)

#Figure out what bacteria are different between 3-day recovery + PEG + FMT or PBS groups (RM vs FRM)
FRMvRM_post_gavage <- genus_pairwise_stools_5plusdpi %>% 
  filter(group1 %in% c( "FRM", "RM") & group2 %in% c("FRM", "RM")) %>% 
  filter(p.adj < .05) %>%
  arrange(p.adj)  %>% 
  distinct(genus) %>% 
  filter(genus != "Unclassified") %>%  #Remove unclassified genus since it's not informative
  pull(genus)
#Create dataframe of RM & FRM mice
frm_rm_subset <-  agg_genus_data_subset_hm %>% 
  filter(group %in% c("FRM", "RM"))
#Create line plot of these genera over time for FRM & RM groups
facet_labels <- FRMvRM_post_gavage
names(facet_labels) <- FRMvRM_post_gavage
line_plot_days <- c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15)
lp_fmt <- line_plot_genus(frm_rm_subset, FRMvRM_post_gavage, line_plot_days, "solid")+
#  scale_x_discrete(limits = c(3, 4, 5:10, 15), breaks = c(3, 4, 5:10, 15), labels = c(3, 4, 5:10, 15)) #Rewrite over -1:10 scale
  scale_x_continuous(limits = c(-1,15), breaks = c(-1:10, 15), labels = c(-1:10, 15)) #Rewrite over -1:10 scale

save_plot(filename = "results/figures/post_CDI_PEG_genus_lineplot_fmt.png", lp_fmt, base_height = 5, base_width = 8)
