source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match post_CDI_PEG Plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8", "7f5f1e") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "CWM", "FRM", "RM", "FMT")
color_labels <- c( "Clind.", "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350", "FMT")

metadata <- metadata %>%
  mutate(day = as.integer(day))  #Day variable (transformed to integer to get rid of decimals on PCoA animation

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Alpha Diversity Analysis-------
# Pull in diversity for alpha diversity analysis using post CDI PEG subset from defined in utilties.R
#Diversity data for all days
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
count_subset(diversity_stools) %>% View()
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

#Plot Shannon over time days -1 to 30 for post CDI PEG subset
shannon_post_cdi_peg_overtime_full <- diversity_data_subset %>%
  filter(group != "FMT") %>% #drop FMTs
  plot_shannon_overtime() +
  scale_x_continuous(breaks = c(-1:10, 15, 20, 25, 30),
                     limits = c(-2,35), #removes day -15 here
                     minor_breaks = c(-2.5:7.5)) +
  scale_y_continuous(limits = c(0,4))+
  labs(x = "Day",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none") #Removing legend to save separately
save_plot("results/figures/post_CDI_PEG_shannon_overtime.png", shannon_post_cdi_peg_overtime_full) #Save full Shannon over time plot without legend


#Plot Shannon over time for first 10 days for post CDI PEG subset
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
  labs(x = "Days Post-Infection",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none") #Removing legend to save separately
save_plot("results/figures/post_CDI_PEG_shannon_stool.png", shannon_post_cdi_peg_overtime_stool) #Save full Shannon over time plot without legend

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

#Plot Sobs overtime 
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
                     limits = c(-15,35)) +
  theme_classic()+
  theme(legend.position = "none") + #Remove legend
  labs(x = "Day",
       y = "Number of Observed Species")
save_plot("results/figures/post_CDI_PEG_sobs_overtime.png", sobs_post_CDI_PEG)

#Sobs oVertime 10 Day Version
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
                     limits = c(-1.5,11)) +
  theme_classic()+
  theme(legend.position = "none") + #Remove legend
  labs(x = "Day",
       y = "Number of Observed Species")
save_plot("results/figures/post_CDI_PEG_sobs_overtime_10d.png", sobs_post_CDI_PEG_10d)

#Sobs over time for stools
x_annotation <- sig_richness_days_stools
y_position <- max(diversity_stools$sobs)+5
label <- kw_label(kw_richness_stools)
sobs_post_CDI_PEG_stool <- diversity_stools %>%
  filter(group != "FMT") %>% #Remove FMTs
  plot_richness_overtime()+
  scale_x_continuous(breaks = c(-1:10, 15, 20, 25, 30),
                     limits = c(-2, 31),
                     minor_breaks = c(-1.5:10.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5)) +
  scale_y_continuous(limits = c(0,122))+
  theme_classic()+
  theme(legend.position = "none") #Remove legend
save_plot("results/figures/post_CDI_PEG_sobs_overtime_stool.png", sobs_post_CDI_PEG_stool)


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
  geom_point()+ theme_classic()
group_legend <- get_legend(group_legend)
save_plot("results/figures/post_CDI_PEG_pcoa_legend.png", group_legend, base_height = 1.2, base_width = 3.2)

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

pcoa_subset_plot_stool <- plot_pcoa(pcoa_post_cdi_peg_stool)+
  labs(x = paste("PCoA 1 (", axis1_stool, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2_stool,"%)", sep = ""))
save_plot(filename = paste0("results/figures/post_CDI_PEG_stool_pcoa.png"), pcoa_subset_plot_stool, base_height = 5, base_width = 5)

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
  filter(sample_type =="stool") %>% #Exclude the other sample types and just perform test on the stools
  mutate(day = fct_relevel(day, "-15", "-1", "0", "1", "2", "3", "4", "5", "6", "7", 
                           "8", "9", "10", "15", "20", "25", "30"))

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


#Plots of the top 20 OTUs that varied across sources at each timepoint----
#Day 1: top 20 OTUs that vary across sources
D1top20_otus <- plot_otus_dx(`sig_otu_day1`[1:20], 1) +#Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D1top20_otus.png", D1top20_otus, base_height = 9, base_width = 7)
#Day 3: top 20 OTUs that vary across sources
D3top20_otus <- plot_otus_dx(`sig_otu_day3`[1:20], 3) + #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D3top20_otus.png", D3top20_otus, base_height = 9, base_width = 7)
#Day 6: top 20 OTUs that vary across sources
D6top20_otus <- plot_otus_dx(`sig_otu_day6`[1:20], 6)+ #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D6top20_otus.png", D6top20_otus, base_height = 9, base_width = 7)
#Day 8: top 20 OTUs that vary across sources
D8top20_otus <- plot_otus_dx(`sig_otu_day8`[1:20], 8)+ #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D8top20_otus.png", D8top20_otus, base_height = 9, base_width = 7)
#Day 10: top 20 OTUs that vary across sources
D10top20_otus <- plot_otus_dx(`sig_otu_day10`[1:20], 10)+ #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D10top20_otus.png", D10top20_otus, base_height = 9, base_width = 7)


#Perform pairwise comparisons for day -15 and -1
# Perform pairwise Wilcoxan rank sum tests for otus that were significantly different across sources of mice on a series of days----
pairwise_day_otu <- function(timepoint, sig_otu_dayX){
  otu_stats <- post_cdi_PEG_subset(agg_otu_data) %>% 
    filter(day == timepoint) %>%
    select(group, otu, agg_rel_abund) %>% 
    group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$group)) %>% tidy())) %>% 
    mutate(median = map(data, get_rel_abund_median_group)) %>% 
    unnest(c(model, median)) %>% 
    ungroup()
  pairwise_stats <- otu_stats %>% 
    filter(otu %in% sig_otu_dayX) %>% 
    group_by(otu) %>% 
    mutate(model=map(data, ~pairwise.wilcox.test(x=.x$agg_rel_abund, g=as.factor(.x$group), p.adjust.method="BH") %>% 
                       tidy() %>% 
                       mutate(compare=paste(group1, group2, sep="-")) %>% 
                       select(-group1, -group2) %>% 
                       pivot_wider(names_from=compare, values_from=p.value)
    )
    ) %>% 
    unnest(model) %>% 
    select(-data, -parameter, -statistic, -p.value) %>% #Get rid of p.value since it's the unadjusted version
    write_tsv(path = paste0("data/process/post_CDI_PEG_otu_stats_day_", timepoint, "_sig.tsv"))
  #Format pairwise stats to use with ggpubr package
  plot_format_stats <- pairwise_stats %>% 
    select(-method, -C, -CWM, -FRM, -RM) %>% 
    group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
    lapply(tidy_pairwise_otu) %>% 
    bind_rows()
  return(plot_format_stats)  
}
##Pairwise test of significant days (1,3,6,8,10) and their corresponding OTUs for stool samples 

otu_day1_stats <- pairwise_day_otu(1, `sig_otu_day1`)
otu_day3_stats <- pairwise_day_otu(3, `sig_otu_day3`)
otu_day6_stats <- pairwise_day_otu(6, `sig_otu_day6`)
otu_day8_stats <- pairwise_day_otu(8, `sig_otu_day8`)
otu_day10_stats <- pairwise_day_otu(10, `sig_otu_day10`)
#Combine pairwise dataframes for all days
otu_pairwise_stools <- rbind(otu_day1_stats, otu_day3_stats, otu_day6_stats, otu_day8_stats, otu_day10_stats)

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
facet_labels <- c("Clind.", "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350") #Create descriptive labels for facets
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
  filter(sample_type =="stool") %>% #Exclude the other sample types and just perform test on the stools
  mutate(day = fct_relevel(day, "-15", "-1", "0", "1", "2", "3", "4", "5", "6", "7", 
                           "8", "9", "10", "15", "20", "25", "30"))

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

#Plots of the top genera that varied across sources at each timepoint----Days 1, 3, 6, 8, 10
D1top_genera <- plot_genus_dx(agg_genus_data_subset, `sig_genus_day1`[1:20], 1) +#Pick top significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D1top_genera.png", D1top_genera, base_height = 9, base_width = 7)
D3top_genera <- plot_genus_dx(agg_genus_data_subset, `sig_genus_day3`[1:20], 1) +#Pick top significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D3top_genera.png", D3top_genera, base_height = 9, base_width = 7)
D6top_genera <- plot_genus_dx(agg_genus_data_subset, `sig_genus_day6`[1:20], 1) +#Pick top significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D6top_genera.png", D6top_genera, base_height = 9, base_width = 7)
D8top_genera <- plot_genus_dx(agg_genus_data_subset, `sig_genus_day8`[1:20], 1) +#Pick top significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D8top_genera.png", D8top_genera, base_height = 9, base_width = 7)
D10top_genera <- plot_genus_dx(agg_genus_data_subset, `sig_genus_day10`[1:20], 1) +#Pick top significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/post_CDI_PEG_D10top_genera.png", D10top_genera, base_height = 9, base_width = 7)

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
facet_labels <- c("Clind.", "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350") #Create descriptive labels for facets
names(facet_labels) <- c("C", "CWM", "FRM", "RM") #values that correspond to group, which is the variable we're faceting by
hm_stool <- hm_plot_genus(agg_genus_data_subset, hm_sig_genera_p_adj, hm_stool_days)+
  scale_x_discrete(breaks = c(-1:10, 15, 20, 25, 30), labels = c(-1:10, 15, 20, 25, 30)) 
save_plot(filename = "results/figures/post_CDI_PEG_genus_heatmap_stools.png", hm_stool, base_height = 14, base_width = 15)

#Plot heat map of genera that were significant in the 5 day subset (repeat with genera?)

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

