source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match 5_days_PEG plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMR")
color_labels <- c( "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 + 10-day recovery")

metadata <- metadata %>%
  mutate(day = as.integer(day))  #Day variable (transformed to integer to get rid of decimals on PCoA animation

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Alpha diversity analysis----

#Subset diversity data to just the 5-day PEG subset:
diversity_data <- five_day_PEG_subset(diversity_data)

group_day_summary <- diversity_data %>%
  group_by(group) %>%
  count(day)

#Plot of shannon diversity at days 1, 4, and 10 when we have sequencing data for 3 groups
shannon_select_days <- diversity_data %>%
  filter(day %in% c(-15, -10, -5, -1, 0, 1, 2, 3, 4, 5, 6, 10, 15, 20, 25, 30)) %>%
  group_by(group, day) %>%
  mutate(median_shannon = median(shannon)) %>% #create a column of median values for each group
  ungroup() %>%
  ggplot(aes(x=group, y =shannon, colour= group))+
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
  geom_jitter(size=2, alpha=0.6, show.legend = FALSE) +
  labs(title=NULL,
       x=NULL,
       y="Shannon Diversity Index")+
  ylim(0, 4)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank()) #Remove x axis ticks
save_plot("results/figures/5_days_PEG_shannon.png", shannon_select_days) #Use save_plot instead of ggsave because it works better with cowplot

#Plot of WMR group over the days we have sequencing data for 3 groups
shannon_WMR <- diversity_data %>%
  filter(group == "WMR") %>%
  group_by(day) %>%
  mutate(median_shannon = median(shannon)) %>% #create a column of median values for each group
  ungroup() %>%
  ggplot(aes(x=group, y =shannon, colour= group))+
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
  geom_jitter(size=2, alpha=0.6, show.legend = FALSE) +
  labs(title=NULL,
       x=NULL,
       y="Shannon Diversity Index")+
  ylim(0, 4)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank()) #Remove x axis ticks
save_plot("results/figures/5_days_PEG_shannon_WMR.png", shannon_WMR) #Use save_plot instead of ggsave because it works better with cowplot

#Plot of sobs (richness) at days 1, 4, and 10 when we have sequencing data for 3 groups
sobs_select_days <- diversity_data %>%
  filter(day %in% c(-15, -10, -5, -1, 0, 1, 2, 3, 4, 5, 6, 10, 15, 20, 25, 30)) %>%
  group_by(group, day) %>%
  mutate(median_sobs = median(sobs)) %>% #create a column of median values for each group
  ungroup() %>%
  ggplot(aes(x=group, y =sobs, colour= group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  #  scale_shape_manual(name=NULL,
  #                     values=shape_scheme,
  #                     breaks=shape_experiment,
  #                     labels=shape_experiment) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median_sobs, ymin = median_sobs), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(size=2, alpha=0.6, show.legend = FALSE) +
  labs(title=NULL,
       x=NULL,
       y="Number of Observed OTUs")+
  ylim(0, 150)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank()) #Remove x axis ticks
save_plot("results/figures/5_days_PEG_richness.png", sobs_select_days) #Use save_plot instead of ggsave because it works better with cowplot


#Plot of sobs (richness) for the WMR group over the days we have sequencing data for 3 groups
sobs_WMR <- diversity_data %>%
  filter(group == "WMR") %>%
  group_by(day) %>%
  mutate(median_sobs = median(sobs)) %>% #create a column of median values for each group
  ungroup() %>%
  ggplot(aes(x=group, y =sobs, colour= group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  #  scale_shape_manual(name=NULL,
  #                     values=shape_scheme,
  #                     breaks=shape_experiment,
  #                     labels=shape_experiment) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median_sobs, ymin = median_sobs), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(size=2, alpha=0.6, show.legend = FALSE) +
  labs(title=NULL,
       x=NULL,
       y="Number of Observed OTUs")+
  ylim(0, 150)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank()) #Remove x axis ticks
save_plot("results/figures/5_days_PEG_richness_WMR.png", sobs_WMR) #Use save_plot instead of ggsave because it works better with cowplot


#Distance matrix of 5_day_PEG PCoA subset----
dist <- read_dist("data/process/5_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")

#Plot PCoA data----
#PCoA plot that combines the 2 experiments and save the plot----

#Read in pcoa loadings and axes for 5_day_PEG PCoA subset
#Pull 5_Day_PEG subset of PCoA data
pcoa_5_day_PEG <- read_tsv("data/process/5_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") 
select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  right_join(diversity_data, by= "unique_label") %>% #merge metadata and PCoA data frames (This drops some of our 16S data for early timepoints)
  mutate(day = as.integer(day)) %>% #Day variable (transformed to integer to get rid of decimals on PCoA animation
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

#Pull axes from loadings file
pcoa_axes_5_day_PEG <- read_tsv("data/process/5_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot <- plot_pcoa(pcoa_5_day_PEG)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))


save_plot(filename = paste0("results/figures/5_Day_PEG_PCoA.png"), pcoa_subset_plot, base_height = 5, base_width = 4.5)


#Remove legend
pcoa_plot_time <- plot_pcoa(pcoa_5_day_PEG)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))+
  theme(legend.position = "none")+ #remove legend
  facet_wrap(~ day)

#PCoAs of select timepoints of interst

#Animation of PCoA plot over time for all sequenced samples ----
#Source: Will Close's Code Club from 4/12/2020 on plot animation
pcoa_animated <- plot_pcoa(pcoa_data)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))+
  labs(title = 'Day: {frame_time}') + #Adds time variable to title
  transition_time(day)+  #Day variable used to cycle through time on animation
  shadow_mark() #Shows previous timepoints

# Implement better frames per second for animation
pcoa_gif <- animate(pcoa_animated, duration = 10, fps = 10,
                    res = 150, width = 20, height = 20, unit = "cm")

# Save as gif file
anim_save(animation = pcoa_gif, filename = 'results/5_days_PEG_pcoa_over_time.gif')


#OTU analysis----
#11/4/20 Note this was implemented for only plates1_2 of 16S sequenced samples.
#Need to update to include all timepoints/tissues now that we have all the sequence data

#Figure out which days we have sequencing data for from the 3 groups:
test <- pcoa_data %>% group_by(group) %>% count(day)
#Days that we have data for all 3 groups: 1, 4, 10
test_days <- c(1, 4, 10)

#Function to test for differences across groups at the OTU level for specific timepoints
#Function to test at the otu level:
kruskal_wallis_otu <- function(timepoint){
  otu_stats <- agg_otu_data %>%
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
    select(otu, statistic, p.value, parameter, method, C, WM, WMR) %>%
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
    arrange(p.value.adj) %>%
    write_tsv(path = paste0("data/process/5_days_PEG_otu_stats_day_", timepoint, ".tsv"))
}

# Perform kruskal wallis tests at the otu level for all days of the experiment that were sequenced----
for (d in test_days){
  kruskal_wallis_otu(d)
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/5_days_PEG_otu_stats_day_", d, ".tsv"))
  name <- paste("sig_otu_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, otu))
}

#OTUs that varied across treatment groups and were shared across days 1, 4 and 10
shared_sig_otus_d1_4_10 <- intersect_all(sig_otu_day1, sig_otu_day4, sig_otu_day10)
#2 OTUs
shared_sig_otus_d1_4 <- intersect_all(sig_otu_day1, sig_otu_day4)
#15 OTUs
shared_sig_otus_d1_10 <- intersect_all(sig_otu_day1, sig_otu_day10)
#2 OTUs
shared_sig_otus_d4_10 <- intersect_all(sig_otu_day4, sig_otu_day10)
#39 OTUs

#Only Lachnospiraceae (OTU 134) and Lachnospiraceae (OTU 126) varies across all 3 timepoints
#25 OTUs vary on day 1: 2 OTUs shared with day 10
#123 OTUs on day 4: 16 OTUs shared with day 1
#47 OTUs on day 10: 23/25 shared with day 4

#Function to plot a list of OTUs across sources of mice at a specific timepoint:
#Arguments: otus = list of otus to plot; timepoint = day of the experiment to plot
plot_otus_dx <- function(otus, timepoint){
  agg_otu_data %>%
    filter(otu %in% otus) %>%
    filter(day == timepoint) %>%
    mutate(agg_rel_abund = agg_rel_abund + 1/4000) %>% # 4,000 is 2 times the subsampling parameter of 2000
    ggplot(aes(x= otu_name, y=agg_rel_abund, color=group))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    geom_hline(yintercept=1/2000, color="gray")+
    stat_summary(fun = 'median',
                 fun.max = function(x) quantile(x, 0.75),
                 fun.min = function(x) quantile(x, 0.25),
                 position = position_dodge(width = 1)) +
    labs(title=NULL,
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "none",
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Plots of the relative abundances of OTUs that significantly varied across sources of mice from day -1 to day 1----
otus_d1 <- plot_otus_dx(sig_otu_day1, 1)+
  ggtitle("Day 1 post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_otus_d1.png", otus_d1, base_height = 7, base_width = 8)

otus_d4 <- plot_otus_dx(sig_otu_day4, 4)+
  ggtitle("Day 4 post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_otus_d4.png", otus_d4, base_height = 7, base_width = 8)

otus_d10 <- plot_otus_dx(sig_otu_day10, 10)+
  ggtitle("Day 10 post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_otus_d10.png", otus_d10, base_height = 7, base_width = 8)

#Examine impacts of clindamycin and PEG3350 treatments on bacterial OTUs----

#Examine changes that happen after clindamycin treatment (baseline day -5 versus day 1)
C_dn5_d1_pairs <- agg_otu_data %>%
  filter(group == "C" & otu == "Bacteroides (OTU 1)") %>% #Limit to group "C" and randomly pick an OTU just to figure out what mice have sequence data
  filter(day == -5 | day == 1) %>%
  filter(duplicated(unique_mouse_id)) %>% #Pull mouse ids with sequence data for both day -1 and day 0
  pull(unique_mouse_id) #6 mice

#Dataframe for statistical test at the OTU level
C_paired_otu <- agg_otu_data %>%
  filter(unique_mouse_id %in% C_dn5_d1_pairs) %>% #Only select pairs with data for day -1 & day 0
  filter(day == -5 | day == 1) %>% #Experiment days that represent initial community and community post clindamycin treatment
  mutate(day = as.factor(day)) %>%
  select(day, otu, agg_rel_abund)

#Wilcoxon signed rank test for all day -1, day 0 pairs at the OTU level:
otus_C_pairs <- C_paired_otu %>%
  group_by(otu) %>%
  nest() %>%
  mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>%
  mutate(median = map(data, get_rel_abund_median_day)) %>%
  unnest(c(model, median)) %>%
  ungroup()

#Adjust p-values for testing multiple OTUs
otus_C_pairs_stats_adjust <- otus_C_pairs %>%
  select(otu, statistic, p.value, method, alternative, `-5`, `1`) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv(path = "data/process/5_days_PEG_otu_Cgroup_dn5to0.tsv")

#Make a list of significant OTUs impacted by clindamycin treatment----
C_sig_otu_pairs <- pull_significant_taxa(otus_C_pairs_stats_adjust, otu)
# 0 OTUs
C_sig_otu_pairs_top10 <- C_sig_otu_pairs[1:10]

C_top_OTUs <- otus_C_pairs_stats_adjust %>%
  arrange(p.value)
C_top10_OTUs <- head(C_top_OTUs, 10) %>% pull(otu)

#Examine OTUs that change after 5-day PEG3350 treatment (baseline day -5 versus day 1)
WM_dn5_d1_pairs <- agg_otu_data %>%
  filter(group == "WM" & otu == "Bacteroides (OTU 1)") %>% #Limit to group "C" and randomly pick an OTU just to figure out what mice have sequence data
  filter(day == -5 | day == 1) %>%
  filter(duplicated(unique_mouse_id)) %>% #Pull mouse ids with sequence data for both day -1 and day 0
  pull(unique_mouse_id) #9 mice

#Dataframe for statistical test at the OTU level
WM_paired_otu <- agg_otu_data %>%
  filter(unique_mouse_id %in% WM_dn5_d1_pairs) %>% #Only select pairs with data for day -1 & day 0
  filter(day == -5 | day == 1) %>% #Experiment days that represent initial community and community post clindamycin treatment
  mutate(day = as.factor(day)) %>%
  select(day, otu, agg_rel_abund)

#Wilcoxon signed rank test for all day -1, day 0 pairs at the OTU level:
otus_WM_pairs <- WM_paired_otu %>%
  group_by(otu) %>%
  nest() %>%
  mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>%
  mutate(median = map(data, get_rel_abund_median_day)) %>%
  unnest(c(model, median)) %>%
  ungroup()

#Adjust p-values for testing multiple OTUs
otus_WM_pairs_stats_adjust <- otus_WM_pairs %>%
  select(otu, statistic, p.value, method, alternative, `-5`, `1`) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv(path = "data/process/5_days_PEG_otu_WMgroup_dn5to0.tsv")

#Make a list of significant OTUs impacted by PEG3350 treatment----
WM_sig_otu_pairs <- pull_significant_taxa(otus_WM_pairs_stats_adjust, otu)
# 0 OTUs
WM_sig_otu_pairs_top10 <- WM_sig_otu_pairs[1:10]

WM_top_OTUs <- otus_WM_pairs_stats_adjust %>%
  arrange(p.value)
WM_top10_OTUs <- head(WM_top_OTUs, 10) %>% pull(otu)

#Compare OTUs impacted by clindamycin and PEG3350
WM_C_otus <- intersect_all(C_top10_OTUs, WM_top10_OTUs)
#4 OTUs overlap: Enterobacteriaceae (OTU 3) and Porphyromonadaceae (OTUs 8, 11, 16)
