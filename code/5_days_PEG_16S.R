source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match 5_days_PEG plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMR")
color_labels <- c( "Clind.", "5-day PEG", "5-day PEG + Clind.", "5-day PEG + 10-day recovery")
#Need to create an additional color scheme with 6 colors (or consider keeping colors and doing open/closed for mock challenged mice)
#See 5_days_PEG_histology_scores.R for how mock challenged mice were presented
#Define shape scheme based on Infected status----
shape_scheme <- c(1, 19)
shape_infected <- c("no", "yes")

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Alpha diversity analysis----

#Subset diversity data to just the 5-day PEG subset:
diversity_subset <- five_day_PEG_subset(diversity_data)
#five_day_PEG_subset() will exclude mock challenged mice (group = WMN or CN)

#Create subset dataframes of the 5-days PEG diversity data for just stool samples, tissues. 
diversity_stools <- subset_stool(diversity_subset) %>%
  mutate(day = as.numeric(day))
diversity_tissues <- subset_tissue(diversity_subset) %>%
  mutate(day = as.numeric(day))

#Figure out how many samples we have per group per day for each subset
num_stool <- count_subset(diversity_stools) #Number of stool samples per group per day
num_tissue <- count_subset(diversity_tissues) #Number of tissue samples per group per day

#Experimental days to analyze with the Kruskal-Wallis test (timepoints with 16S data for at least 3 groups)
#Baseline (before treatment) for WMR is day -15. For C, WM, and WMC baseline is day -5
stool_test_days <- c(-5, -1, 0, 1, 2, 3, 4, 5, 6, 10, 30)
tissue_test_days <- c(6, 30) #Only 2 days with samples from at least 3 groups

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
    select(-data) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/5_days_PEG_shannon_stats_", subset_name, "_subset.tsv"))
}
#Test with shannon for stool subset
kw_shannon_stools <- kruskal_wallis_shannon(diversity_stools, stool_test_days, "stools")
sig_shannon_days_stools <- pull_sig_days(kw_shannon_stools)
#Test with shannon for tissue subset
kw_shannon_tissues <- kruskal_wallis_shannon(diversity_tissues, tissue_test_days, "tissues")
sig_shannon_days_tissues <- pull_sig_days(kw_shannon_tissues)

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
    select(-data) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/5_days_PEG_richness_stats_", subset_name, "_subset.tsv"))
}
#Test with richness for stool subset
kw_richness_stools <- kruskal_wallis_richness(diversity_stools, stool_test_days, "stools")
sig_richness_days_stools <- pull_sig_days(kw_richness_stools)
#Test with richness for tissue subset
kw_richness_tissues <- kruskal_wallis_richness(diversity_tissues, tissue_test_days, "tissues")
sig_richness_days_tissues <- pull_sig_days(kw_richness_tissues)

#Shannon and richness plots for stool samples from the 5-days PEG subset----
#Plot Shannon diversity over time for the subset of stool samples (excluding mock challenged mice):
#Statistical annotation labels:
x_annotation <- sig_shannon_days_stools
y_position <- max(diversity_stools$shannon)+0.15
label <- kw_label(kw_shannon_stools)
#Plot
shannon_stools <- plot_shannon_overtime(diversity_stools) +
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_shannon_stools.png", shannon_stools, base_height = 4, base_width = 9, base_aspect_ratio = 2)

#Plot richness over time for the subset of stool samples
#Statistical annotation labels:
x_annotation <- sig_richness_days_stools
y_position <- max(diversity_stools$sobs)+5
label <- kw_label(kw_richness_stools)
#Plot
richness_stools <- plot_richness_overtime(diversity_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_richness_stools.png", richness_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Shannon and richness plots for stool samples from the 5-days PEG subset----
#Statistical annotation labels:
x_annotation <- sig_shannon_days_tissues
y_position <- max(diversity_tissues$shannon)+ 0.05
label <- kw_label(kw_shannon_tissues)
#Plot
shannon_tissues <- plot_shannon_overtime(diversity_tissues) +
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_shannon_tissues.png", shannon_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Statistical annotation labels:
x_annotation <- sig_richness_days_tissues
y_position <- max(diversity_tissues$sobs)+5
label <- kw_label(kw_richness_tissues)
#Plot
richness_tissues <- plot_richness_overtime(diversity_tissues) +
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_richness_tissues.png", richness_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Distance matrix of 5_day_PEG PCoA subset----
dist <- read_dist("data/process/5_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")

#Plot PCoA data----

#Read in pcoa loadings and axes for 5_day_PEG PCoA subsets (stools and tissues). Excludes mock samples

#Pull 5_Day_PEG subset PCoA of stool samples
pcoa_5_day_PEG_stool <- read_tsv("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  mutate(day = as.integer(day)) %>% #Transform day into continuous
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

pcoa_axes_5_day_PEG_stool <- read_tsv("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot_stool <- plot_pcoa(pcoa_5_day_PEG_stool)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  theme(legend.position = "none")
save_plot(filename = paste0("results/figures/5_days_PEG_stool_PCoA.png"), pcoa_subset_plot_stool, base_height = 5, base_width = 5)

##Animation of PCoA plot: Stool Subset--
pcoa_animated_stool <- plot_pcoa(pcoa_5_day_PEG_stool)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  labs(title = 'Day: {frame_time}') + #Adds time variable to title
  transition_time(day)+  #Day variable used to cycle through time on animation
  shadow_mark() #Shows previous timepoints

# Implement better frames per second for animation
pcoa_gif_stool <- animate(pcoa_animated_stool, duration = 6, fps = 10,
                          res = 150, width = 20, height = 20, unit = "cm")

# Save as gif file
anim_save(animation = pcoa_gif_stool, filename = 'results/5_days_PEG_pcoa_over_time_stools.gif')

#Create stand alone legend
group_legend <- pcoa_5_day_PEG_stool  %>%
  ggplot(aes(x = axis1, y = axis2, color = group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  geom_point()+ theme_classic()
group_legend <- get_legend(group_legend)
save_plot("results/figures/5_days_PEG_pcoa_legend.png", group_legend, base_height = 1, base_width = 2.3)

#Tissue subset
pcoa_5_day_PEG_tissues <- read_tsv("data/process/5_day_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and use left_join to keep all samples in pcoa data frame
  mutate(day = as.integer(day)) %>% #Transform day into integer
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

#Pull axes from loadings file
pcoa_axes_5_day_PEG_tissues <- read_tsv("data/process/5_day_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG_tissues %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG_tissues %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot_tissue <- plot_pcoa(pcoa_5_day_PEG_tissues)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  scale_alpha_continuous(range = c(.3, 1),
                         breaks= c(4, 6, 20, 30),
                         labels=c(4, 6, 20, 30))+
  theme(legend.position = "none")
save_plot(filename = paste0("results/figures/5_days_PEG_tissues_PCoA.png"), pcoa_subset_plot_tissue, base_height = 5, base_width = 5)

#PCoA faceted over time
pcoa_plot_time <- plot_pcoa(pcoa_5_day_PEG_tissues)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))+
  theme(legend.position = "none")+ #remove legend
  facet_wrap(~ day)

##Animation of PCoA plot: Tissue Subset--
pcoa_animated_tissues <- plot_pcoa(pcoa_5_day_PEG_tissues)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  labs(title = 'Day: {frame_time}') + #Adds time variable to title
  transition_time(day)+  #Day variable used to cycle through time on animation
  shadow_mark() #Shows previous timepoints

# Implement better frames per second for animation
pcoa_gif_tissue <- animate(pcoa_animated_tissues, duration = 6, fps = 10,
                    res = 150, width = 20, height = 20, unit = "cm")

# Save as gif file
anim_save(animation = pcoa_gif_tissue, filename = 'results/5_days_PEG_pcoa_over_time_tissues.gif')

#Genus level analysis----
#Subset genus data to just the 5-day PEG subset and separate by sample type (stools versus tissues)
genus_subset <- five_day_PEG_subset(agg_genus_data)
#Create subset dataframes of the 5-days PEG diversity data for just stool samples, tissues. 
genus_stools <- subset_stool(genus_subset)
genus_tissues <- subset_tissue(genus_subset)
#The above subsets exclude mock challenged mice (group = WMN or CN)
#Also create dataframes of diversity data that includes mock challenged mice (WMN and C), separated into stool and tissue samples
genus_mock_stools <- subset_stool(add_mocks(genus_subset, agg_genus_data))
genus_mock_tissues <- subset_tissue(add_mocks(genus_subset, agg_genus_data))

#Function to test for differences across groups at the genus level for specific timepoints
#Function to test at the genus level:
#Arguments:
# timepoint = day of the experiment
#sample_df = subset dataframe of just stool or tissue samples
#sample_type = "stool" or "tissue" to be included in filename
kruskal_wallis_genus <- function(timepoint, sample_df, sample_type){
  genus_stats <- sample_df %>%
    filter(day == timepoint) %>%
    select(group, genus, agg_rel_abund) %>%
    group_by(genus) %>%
    nest() %>%
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$group)) %>% tidy())) %>%
    mutate(median = map(data, get_rel_abund_median_group)) %>%
    unnest(c(model, median)) %>%
    ungroup()
  #Adjust p-values for testing multiple Genera
  genus_stats_adjust <- genus_stats %>%
    select(-data) %>% #Keep everything but the data column
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
    arrange(p.value.adj) %>%
    write_tsv(path = paste0("data/process/5_days_PEG_genus_stats_day_", timepoint, "_", sample_type, ".tsv"))
}

#Create empty data frame to combine stat dataframes for all days that were tested
kw_genus_stools <- data.frame(genus=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                              WM =double(),C =double(),WMR =double(),WMC=double(),
                              p.value.adj=double(),day=double())

# Perform kruskal wallis tests at the genus level for the stool samples----
for (d in stool_test_days){
  kruskal_wallis_genus(d, genus_stools, "stools")
  #Make a list of significant genus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/5_days_PEG_genus_stats_day_", d, "_stools.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_genus_stools_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, genus))
  kw_genus_stools <- add_row(kw_genus_stools, stats)  #combine all the dataframes together
}

#Save combined Kruskal-Wallis tests at the genus level
kw_genus_stools %>% write_tsv("data/process/5_days_PEG_genus_group_stools.tsv")

#Create empty data frame to combine stat dataframes for all days that were tested
kw_genus_tissues <- data.frame(genus=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                               WM =double(),C =double(),WMR =double(),WMC=double(),
                               p.value.adj=double(),day=double())
# Perform kruskal wallis tests at the genus level for the tissue samples----
for (d in tissue_test_days){
  kruskal_wallis_genus(d, genus_tissues, "tissues")
  #Make a list of significant genus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/5_days_PEG_genus_stats_day_", d, "_tissues.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_genus_tissues_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, genus))
  kw_genus_tissues <- add_row(kw_genus_tissues, stats)  #combine all the dataframes together
}

#Save combined Kruskal-Wallis tests at the genus level
kw_genus_tissues %>%  write_tsv("data/process/5_days_PEG_genus_group_tissues.tsv")

#Examine genera that varied between treatment groups across multiple days
#Stool samples genera
stool_test_days #Check days tested. Look at all timepoints post challenge
shared_sig_stools_genus_d1tod30 <- intersect_all(sig_genus_stools_day1, sig_genus_stools_day2,                
                                                sig_genus_stools_day3, sig_genus_stools_day4, 
                                                sig_genus_stools_day5, sig_genus_tissues_day6,
                                                sig_genus_stools_day10, sig_genus_stools_day30) #fill in different days to compare
shared_sig_stools_genus_d1tod30 #only "Peptostreptococcaceae Unclassified"
#D1 to 10 post challenge
shared_sig_stools_genus_d1tod10 <- intersect_all(sig_genus_stools_day1, sig_genus_stools_day2,                
                                                 sig_genus_stools_day3, sig_genus_stools_day4, 
                                                 sig_genus_stools_day5, sig_genus_tissues_day6,
                                                 sig_genus_stools_day10)
shared_sig_stools_genus_d1tod10 #only "Peptostreptococcaceae Unclassified"
#D1 to 6 post challenge
shared_sig_stools_genus_d1tod6 <- intersect_all(sig_genus_stools_day1, sig_genus_stools_day2,                
                                                 sig_genus_stools_day3, sig_genus_stools_day4, 
                                                 sig_genus_stools_day5, sig_genus_tissues_day6)
shared_sig_stools_genus_d1tod6 #4 genera

#Tissue samples genera
tissue_test_days #Check days tested. d6 and 30
shared_sig_tissues_genus <- intersect_all(sig_genus_tissues_day6, sig_genus_tissues_day30) 
shared_sig_tissues_genus #7 genera including unclassified 
#Drop unclassified bacteria, not informative
shared_sig_tissues_genus <- shared_sig_tissues_genus[1:6]

#Plots of the relative abundances for no more than the top 10 significant genera that varied over the tested timepoints
#Stools
#day1
genus_d1 <- plot_genus_dx(genus_stools, sig_genus_stools_day1[1:10], 1)+
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Stool 1 dpi")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d1.png", genus_d1, base_height = 7, base_width = 8)
#day2
genus_d2 <- plot_genus_dx(genus_stools, sig_genus_stools_day2[1:10], 2)+
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Stool 2 dpi")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d2.png", genus_d2, base_height = 7, base_width = 8)
#day3
genus_d3 <- plot_genus_dx(genus_stools, sig_genus_stools_day3[1:10], 3)+
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Stool 3 dpi")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d3.png", genus_d3, base_height = 7, base_width = 8)
#day4
genus_d4 <- plot_genus_dx(genus_stools, sig_genus_stools_day4[1:10], 4)+
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Stool 4 dpi")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d4.png", genus_d4, base_height = 7, base_width = 8)
#day5
genus_d5 <- plot_genus_dx(genus_stools, sig_genus_stools_day5[1:10], 5)+
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Stool 5 dpi")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d5.png", genus_d5, base_height = 7, base_width = 8)
#day6
genus_d6 <- plot_genus_dx(genus_stools, sig_genus_stools_day6[1:10], 6)+
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Stool 6 dpi")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d6.png", genus_d6, base_height = 7, base_width = 8)
#day10
genus_d10 <- plot_genus_dx(genus_stools, sig_genus_stools_day10[1:10], 10)+
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Stool 10 dpi")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d10.png", genus_d10, base_height = 7, base_width = 8)
#day30
genus_d30 <- plot_genus_dx(genus_stools, sig_genus_stools_day30[1:10], 30)+
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Stool 30 dpi")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d30.png", genus_d30, base_height = 7, base_width = 8)

#Plots of the relative abundances for no more than the top 10 significant genera that varied over the tested timepoints
#Tissues
#day6
genus_tissues_d6 <- plot_genus_dx(genus_tissues, sig_genus_tissues_day6[1:10], 6)+
  ggtitle("Tissue 6 dpi")+ #Title plot
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_tissues_d6.png", genus_tissues_d6, base_height = 7, base_width = 8)
#day30 
genus_tissues_d30 <- plot_genus_dx(genus_tissues, sig_genus_tissues_day6[1:10], 30)+
  ggtitle("Tissue 30 dpi")+ #Title plot
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_tissues_d30.png", genus_tissues_d30, base_height = 7, base_width = 8)

#List of the 10 significant genera in stool samples that are significant over the most timepoints
top_10_sig_genus <- kw_genus_stools %>% 
  filter(p.value.adj < 0.05) %>% #Select only significant p-values
  group_by(genus) %>% 
  tally() %>% 
  filter(n > 1) %>% #select genera that vary over multiple timepoints = n > 1
  arrange(desc(n)) %>% #Rank by how many time each genera shows up (
  filter(!genus == "Unclassified") %>% #Remove this genus since it's not informative and could contain multiple unclassified genera
  head(10) %>% 
  pull(genus)
#Heatmap of the 10 significant genera in stool samples that are significant over the most timepoints----
hm_stool_days <- diversity_stools %>% distinct(day) %>% pull(day) #Have more timepoints than we tested (some days we only have sequenced samples for 1 group)
facet_labels <- color_labels #Create descriptive labels for facets
names(facet_labels) <- c("C", "WM", "WMC", "WMR") #values that correspond to group, which is the variable we're faceting by
hm_stool <- hm_plot_genus(genus_stools %>% 
                            #Reorder list of genera to match number of days they were significantly different, Bacteroides was the top
                            mutate(genus = fct_relevel(genus, rev(top_10_sig_genus))),
                          top_10_sig_genus, hm_stool_days)+
  scale_x_discrete(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30), labels = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_genus_heatmap_stools.png", hm_stool, base_height = 14, base_width = 15)

#Lineplots of the top 6 significant genera in stool samples----
lp_stool_days <- diversity_stools %>% distinct(day) %>% 
  filter(!day %in% c("-15", "-10", "-5", "-4", "-2", "15", "20", "30")) %>% #Focus on day -1 through 10 timepoints
  pull(day)
facet_labels <- top_10_sig_genus[1:6] #Pick just the top 6
names(facet_labels) <- top_10_sig_genus[1:6] #Pick just the top 6
lp_stool <- line_plot_genus(genus_stools, top_10_sig_genus[1:6], lp_stool_days, "solid")
save_plot(filename = "results/figures/5_days_PEG_genus_lineplot_stools.png", lp_stool, base_height = 5, base_width = 8)

#List of the 10 significant genera in tissue samples that are significant over the most timepoints
sig_genus_tissues<- kw_genus_tissues %>% 
  filter(p.value.adj < 0.05) %>% #Select only significant p-values
  filter(!genus == "Unclassified")  #Remove this genus since it's not informative and could contain multiple unclassified genera
#Find genera that are different between groups on day 6 and 30
sig_genus_tissues_multi_day <- sig_genus_tissues %>%  
  group_by(genus) %>% 
  tally() %>% #Count how many times each genus appears
  filter(n > 1) %>% #select genera that vary over day 6 & day 30
  arrange(desc(n)) %>% #Rank by how many time each genera shows up (multiple timepoints = n > 1)
  head(10) %>% 
  pull(genus) #6 genera
#Find 4 more genera based on lowest adjusted p-values (excluding 6 genera already pulled)
sig_genus_tissue_top_p <- sig_genus_tissues %>% 
  arrange(p.value.adj) %>% 
  filter(!genus %in% sig_genus_tissues_multi_day) %>% 
  head(4) %>% 
  pull(genus)
#Combine 2 lists of significant genera to create final list of 10
top_sig_genus_tissues <- c(sig_genus_tissues_multi_day, sig_genus_tissue_top_p)
#Create heatmap of significant genera for tissue samples----
hm_tissue_days <- diversity_tissues %>% distinct(day) %>% pull(day)
hm_tissues <- hm_plot_genus(genus_tissues%>% 
                              #Reorder list of genera to match number of days they were significantly different, Bacteroides was the top
                              mutate(genus = fct_relevel(genus, rev(top_sig_genus_tissues))), 
                            top_sig_genus_tissues, hm_tissue_days)+
  scale_x_discrete(breaks = c(4, 6, 20, 30), labels = c(4, 6, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_genus_heatmap_tissues.png", hm_tissues, base_height = 14, base_width = 15)

#Lineplots of the top 6 significant genera in tissue samples----
lp_tissue_days <- diversity_tissues %>% distinct(day) %>% 
  filter(day %in% c("4", "6", "20", "30")) %>% #Focus on day -1 through 10 timepoints
  pull(day)
facet_labels <- top_sig_genus_tissues[1:6] #Pick just the top 6
names(facet_labels) <- top_sig_genus_tissues[1:6] #Pick just the top 6
lp_tissue <- line_plot_genus(genus_tissues, top_sig_genus_tissues[1:6], lp_tissue_days, "solid")+
  scale_x_continuous(limits = c(3.5,30.5), breaks = c(4, 6, 20, 30), labels = c(4, 6,20, 30))#Change scale since we have less timepoints for tissues
save_plot(filename = "results/figures/5_days_PEG_genus_lineplot_tissues.png", lp_tissue, base_height = 5, base_width = 8)

#List of genera that overlap between top genera for stool and tissue samples----
top_genera <- intersect_all(top_10_sig_genus, top_sig_genus_tissues)
#5 genera: Bacteroides, Clostridiales Unclassified, Ruminococcaceae Unclassified, 
#Peptostreptococcaceae Unclassified, Acetatifactor

#OTU analysis----
#Subset otu data to just the 5-day PEG subset and separate by sample type (stools versus tissues)
otu_subset <- five_day_PEG_subset(agg_otu_data)
#Create subset dataframes of the 5-days PEG diversity data for just stool samples, tissues. 
otu_stools <- subset_stool(otu_subset)
otu_tissues <- subset_tissue(otu_subset)
#The above subsets exclude mock challenged mice (group = WMN or CN)
#Also create dataframes of diversity data that includes mock challenged mice (WMN and C), separated into stool and tissue samples
otu_mock_stools <- subset_stool(add_mocks(otu_subset, agg_otu_data))
otu_mock_tissues <- subset_tissue(add_mocks(otu_subset, agg_otu_data))

#Transform day variable into integer to use continuous scale on plots of OTUs over time
otu_stools <- otu_stools %>% 
  mutate(day = as.integer(day)) 
otu_tissues <- otu_tissues %>% 
  mutate(day = as.integer(day))

#Examine C. difficile otu over time----
peptostrep_stools <- otu_over_time("Peptostreptococcaceae (OTU 12)", otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_otu_peptostreptococcaceae_stools.png", peptostrep_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
peptostrep_tissues <- otu_over_time("Peptostreptococcaceae (OTU 12)", otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_otu_peptostreptococcaceae_tissues.png", peptostrep_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Bacteroides otu over time----
bacteroides_stools <- otu_over_time("Bacteroides (OTU 1)", otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_otu_bacteroides_stools.png", bacteroides_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
bacteroides_tissues <- otu_over_time("Bacteroides (OTU 1)", otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_otu_bacteroides_tissues.png", bacteroides_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Enterobacteriaceae (OTU 2) over time----
entero2_stools <- otu_over_time("Enterobacteriaceae (OTU 2)", otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_otu_enterobacteriaceae_stools.png", entero2_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
entero2_tissues <- otu_over_time("Enterobacteriaceae (OTU 2)", otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_otu_enterobacteriaceae_tissues.png", entero2_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Porphyromonadaceae otu over time----
porph_stools <- otu_over_time("Porphyromonadaceae (OTU 8)", otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_otu_porphyromonadaceae8_stools.png", porph_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
porph_tissues <- otu_over_time("Porphyromonadaceae (OTU 8)", otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_otu_porphyromonadaceae8_porph_tissues.png", porph_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine changes that happen in WMR group "5-day PEG 3350 + 10-day recovery" 
#post C. difficile challenge (day 1 versus day 8)---- 
#Day 8 because that is when the group median stabilizes
WMR_pairs <- genus_stools %>%
  filter(group == "WMR" & genus == "Bacteroides") %>% #Limit to group and randomly pick a genus just to figure out what mice have sequence data
  filter(day == 1 | day == 8) %>%
  filter(duplicated(unique_mouse_id)) %>% #Pull mouse ids with sequence data for both day -1 and day 0
  pull(unique_mouse_id) #9 mice

#Wilcoxon signed rank test for all day 1, day 8 pairs at the genus level:
genus_WMR_pairs <- genus_stools %>%
  filter(unique_mouse_id %in% WMR_pairs) %>% 
  filter(day == 1 | day == 8) %>% #Select timepoints to test
  group_by(genus) %>%
  nest() %>%
  mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>%
  mutate(median = map(data, get_rel_abund_median_day)) %>%
  unnest(c(model, median)) %>%
  ungroup()

#Adjust p-values for testing multiple genera
genus_WMR_pairs_stats_adjust <- genus_WMR_pairs %>%
  select(genus, statistic, p.value, method, alternative, `1`, `8`) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv(path = "data/process/5_days_PEG_genus_WMR_paired.tsv")

#Create list of genera to plot for WMR group
genus_WMR_pairs_stats_adjust %>% filter(p.value.adj < 0.05) #No p-values survive multiple hypothesis correction
#Look at top 10 genera by unadjusted p-values
WMR_genus <- genus_WMR_pairs_stats_adjust %>% 
  arrange(p.value) %>% 
  filter(p.value < 0.1) %>% #Use 0.1 as p.value cutoff (no bacteria survive FDR correction)
  filter(!genus == "Unclassified") %>% #Remove this genera since it is uninformative
  slice_head(n=9) %>% 
  pull(genus)
#Add C. diff to the list of genera to plot, comment out to leave out, primarily undectable in WMR mice
#WMR_genus <- c(WMR_genus, "Peptostreptococcaceae Unclassified") 

#Create lineplots of 10 genera
lp_stool_days <- diversity_stools %>% distinct(day) %>% 
  filter(!day %in% c("-15", "-10", "-5", "-4", "-2")) %>% #Focus on day -1 through 30 timepoints
  pull(day)
facet_labels <- WMR_genus 
names(facet_labels) <- WMR_genus 
lp_stool_WMR <- genus_stools %>% 
    filter(group == "WMR") %>% #Only plot WMR mice
    filter(genus %in% WMR_genus) %>% #Select only genera of interest
    mutate(genus = fct_relevel(genus, WMR_genus)) %>% #Reorder genera to match order of genera of interest
    filter(day %in% lp_stool_days) %>% #Select only timepoints of interest
    mutate(day = as.integer(day)) %>% 
    group_by(group, genus, day) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_line(aes(x = day, y=median, color=group), linetype = "solid")+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    scale_x_continuous(limits = c(-1.5,31), breaks = c(0, 5, 10, 15, 20, 25, 30), labels = c(0, 5, 10, 15, 20, 25, 30))+
    scale_y_continuous(trans = "log10", limits = c(1/10900, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    geom_hline(yintercept=1/1000, color="gray")+ #Represents limit of detection
    labs(title=NULL,
         x="Days Post-Infection",
         y="Relative abundance (%)")+
    facet_wrap(~genus, nrow = 2, labeller = label_wrap_gen(width = 10))+
    theme_classic()+
    theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
          text = element_text(size = 16),
          strip.text = element_text(face = "italic"),
          plot.title = element_markdown(hjust = 0.5), #Have only the genera names show up as italics
          legend.position = "None")
save_plot(filename = "results/figures/5_days_PEG_genus_lineplot_stools_WMR.png", lp_stool_WMR, base_height = 4, base_width = 8.5)

#Heatmap of 5-day PEG + 10 day recovery (WMR) group----
genus_WMR_stools  <-  genus_stools %>% 
  filter(group == "WMR")
hm_WMR_stool_days <- c("1", "2", "3", "4",
                       "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")
hm_WMR_stool <- genus_WMR_stools %>%
    mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
    mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                             "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the groups  
    filter(genus %in% WMR_genus) %>%
    mutate(genus = fct_relevel(genus, rev(WMR_genus))) %>% #Rearrange order of genera to match significance + Peptostreptococcaceae
    filter(day %in% hm_WMR_stool_days) %>% 
    group_by(group, genus, day) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_tile(aes(x = day, y=genus, fill=median))+
    labs(title=NULL,
         x=NULL,
         y=NULL)+
    scale_fill_distiller(trans = "log10",palette = "YlGnBu", direction = 1, name = "Relative \nAbundance",
                         limits = c(1/10000, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          strip.background = element_blank(), #get rid of box around facet_wrap labels
          axis.text.y = element_text(face = "italic"), #Have genera names in italics
          text = element_text(size = 16))+ # Change font size for entire plot+
  scale_x_discrete(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30), labels = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_genus_heatmap_stools_WMR.png", hm_WMR_stool, base_height = 8, base_width = 8)

#Repeat WMR analysis at the OTU level----
#Wilcoxon signed rank test for all day 1, day 8 pairs at the OTU level:
otu_WMR_pairs <- otu_stools %>%
  filter(unique_mouse_id %in% WMR_pairs) %>% 
  filter(day == 1 | day == 8) %>% #Select timepoints to test
  group_by(otu) %>%
  nest() %>%
  mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>%
  mutate(median = map(data, get_rel_abund_median_day)) %>%
  unnest(c(model, median)) %>%
  ungroup()

#Adjust p-values for testing multiple genera
otu_WMR_pairs_stats_adjust <- otu_WMR_pairs %>%
  select(otu, statistic, p.value, method, alternative, `1`, `8`) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv(path = "data/process/5_days_PEG_otu_WMR_paired.tsv")

#Create list of genera to plot for WMR group
otu_WMR_pairs_stats_adjust %>% filter(p.value.adj < 0.05) #No p-values survive multiple hypothesis correction
#Look at top 10 genera by unadjusted p-values
WMR_otu <- otu_WMR_pairs_stats_adjust %>% 
  arrange(p.value) %>% 
  slice_head(n=10) %>% 
  pull(otu)
#Add C. diff to the list of genera to plot
WMR_otu <- c(WMR_otu, "Peptostreptococcaceae (OTU 12)")

#Heatmap of 5-day PEG + 10 day recovery (WMR) group----
otu_WMR_stools  <-  otu_stools %>% 
  filter(group == "WMR")
hm_WMR_stool_days <- c("1", "2", "3", "4",
                       "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")
hm_WMR_stool <- otu_WMR_stools %>%
  mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
  mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                           "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the groups  
  filter(otu %in% WMR_otu) %>%
  filter(day %in% hm_WMR_stool_days) %>% 
  group_by(group, otu_name, day) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  ggplot()+
  geom_tile(aes(x = day, y=otu_name, fill=median))+
  labs(title=NULL,
       x=NULL,
       y=NULL)+
  scale_fill_distiller(trans = "log10",palette = "YlGnBu", direction = 1, name = "Relative \nAbundance",
                       limits = c(1/10000, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5),
        strip.background = element_blank(), #get rid of box around facet_wrap labels
        axis.text.y = element_markdown(), #Have bacteria name in italics
        text = element_text(size = 16))+ # Change font size for entire plot+
  scale_x_discrete(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30), labels = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_otu_heatmap_stools_WMR.png", hm_WMR_stool, base_height = 8, base_width = 8)

#Examine impacts of PEG3350 treatments on bacterial genera----
#Use genus level stool sample data from baseline and day of/after PEG treatment stopped 
#day of PEG treatment stopping was usually harder to get a sample because the mice still had diarrhea

#For WMR group: Baseline (B) is day -15, post-treatment (P) is -4 (1 day after coming off treatment)
#For WM and WMC: Baseline is day -5, post-treatment is 1 (1 day after coming off treatment)
#For C group: Also use Baseline as day -5 and post-treatment is 1 (don't include in statistical test but include as comparison for plot
#Create dataframe with corresponding days renamed as B for baseline or P for post-treatment, select those timepoints
pairwise_genus_stools <- genus_stools %>% 
  mutate(day = case_when(group == "WMR" & day == "-15" ~ "B",
                         group == "WM" & day == "-5" ~ "B",
                         group == "WMC" & day == "-5" ~ "B",
                         group == "C" & day == "-5" ~ "B",
                         group == "WMR" & day == "-4" ~ "P",
                         group == "WM" & day == "1" ~ "P",
                         group == "WMC" & day == "1" ~ "P",
                         group == "C" & day == "1" ~ "P",
                         TRUE ~ day)) %>% 
  filter(day %in% c("B", "P"))

#Create list of unique mouse IDs of WMR and WM groups to include in statistical test, figure out impact of 5-day PEG
#Exclude C, want to see the impact of 5-day PEG treatment
#Exclude WMC, since these mice also had clindamycin treatment
peg_pairwise_mice <- pairwise_genus_stools %>% 
  filter(genus == "Bacteroides") %>% #Pick one genus just to get list of mice
  filter(group %in% c("WMR", "WM")) %>% 
  filter(duplicated(unique_mouse_id)) %>%
  pull(unique_mouse_id) #25 mice total
  
#Wilcoxon signed rank test for all baseline, post-treatment WM & WMR pairs at the genus level:
peg_wilcoxon <- pairwise_genus_stools %>% 
  filter(unique_mouse_id %in% peg_pairwise_mice) %>% #Select WM & WMR with baseline & post-treatment samples
  group_by(genus) %>% 
  nest() %>%
  mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>%
  mutate(median = map(data, get_rel_abund_median_day)) %>%
  unnest(c(model, median)) %>%
  ungroup()
#Adjust p-values for testing multiple genera
peg_wilcoxon_adjust <- peg_wilcoxon %>% 
  select(genus, statistic, p.value, method, `B`, `P`) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>%
  write_tsv(path = "data/process/5_days_PEG_genus_PEG_paired.tsv")
 
#Pull significant genera from adjusted p-value dataframes
peg_wilcoxon_adjust_sig <- pull_significant_taxa(peg_wilcoxon_adjust, genus)
#18 genera significant

#Create list of genera to plot
peg_wilcoxon_adjust_sig_plot <- peg_wilcoxon_adjust %>% 
  filter(p.value < 0.05) %>% 
  #Exclude Unclassified because it's not informative
  #Exclude Peptostreptococcaceae since we challenged mice with that on day 0 and so will also show up on day 1 (post-treatment for WM, WMC, and WMR mice)
  filter(!genus %in% c("Unclassified", "Peptostreptococcaceae Unclassified")) %>% 
  #Select top 14
  head(14) %>% 
  pull(genus)

#Plot of significant genera  that change between baseline and PEG treatment----
#Include C and WMC group for comparison
facet_labels <- c("Baseline", "Post-treatment")
names(facet_labels) <- c("B", "P")
peg_impacted_genera_plot <- pairwise_genus_stools %>% 
  filter(genus %in% peg_wilcoxon_adjust_sig_plot) %>% 
  filter(day %in% c("B", "P")) %>% #Select baseline and post-treatment timepoints
  mutate(genus=fct_relevel(genus, levels=rev(peg_wilcoxon_adjust_sig_plot))) %>% #Reorder list of genera in order of significance
  mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% 
  ggplot(aes(x= genus, y=agg_rel_abund, color=group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  geom_hline(yintercept=1/1000, color="gray")+
  stat_summary(fun = 'median', 
               fun.max = function(x) quantile(x, 0.75), 
               fun.min = function(x) quantile(x, 0.25),
               position = position_dodge(width = 1)) +  
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/1/10000, 1))+
  coord_flip()+
  theme_classic()+
  geom_vline(xintercept = c((1:18) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  facet_wrap(~day, labeller = labeller(day = facet_labels), scales = "fixed", )+
  theme(plot.title=element_text(hjust=0.5),
        text = element_text(size = 16),# Change font size for entire plot
        axis.text.y = element_markdown(face = "italic"), #Make sure genera names are in italics
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"), #Increase spacing between facets
        legend.position = "none") 
save_plot(filename = paste0("results/figures/5_days_PEG_genera_impacted_by_PEG.png"), peg_impacted_genera_plot, base_height = 9, base_width = 9)

#Alternative barbell plot for bacteria that changed with PEG treatment----
#Format data for barbell plot
peg_impacted_genera_bb <- pairwise_genus_stools %>% 
  filter(genus %in% peg_wilcoxon_adjust_sig_plot) %>% 
  filter(day %in% c("B", "P")) %>% #Select baseline and post-treatment timepoints
  mutate(genus=fct_relevel(genus, levels=rev(peg_wilcoxon_adjust_sig_plot))) %>% #Reorder list of genera in order of significance
  group_by(group, genus, day) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  pivot_wider(names_from = day, values_from = median) %>% 
  ungroup
#Create barbell plot
facet_labels <- color_labels
names(facet_labels) <- color_groups
peg_impacted_genera_bb_plot <- peg_impacted_genera_bb %>% 
  ggplot(aes(x= genus, color=group))+
  geom_segment(aes(y= B, yend = P, xend = genus), size = .5,
               arrow= arrow(length=unit(0.30,"cm"), ends="last"))+
  geom_point(aes(y = B), shape = 1, stroke = 1, size = 3) + 
  geom_point(aes(y = P), shape = 16, size = 3) + 
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  geom_hline(yintercept=1/1000, color="gray")+
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/1/10000, 1))+
  coord_flip()+
  theme_classic()+
  facet_wrap(~group, labeller = labeller(group = facet_labels), scales = "fixed")+
  geom_vline(xintercept = c((1:18) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(plot.title=element_text(hjust=0.5),
        text = element_text(size = 16),# Change font size for entire plot
        axis.text.y = element_markdown(face = "italic"), #Make sure genera names are in italics
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"), #Increase spacing between facets
        legend.position = "none") 
save_plot(filename = paste0("results/figures/5_days_PEG_genera_impacted_by_PEG_barbell.png"), peg_impacted_genera_bb_plot, base_height = 9, base_width = 8)


#Examine genera of interest over time----
genus_stools_t <- genus_stools %>% 
  mutate(day = as.integer(day))  #Day variable transformed to integer for over time plots
genus_tissues_t <- genus_tissues %>% 
  mutate(day = as.integer(day))  #Day variable transformed to integer for over time plots  

#Examine Bacteroides over time
bacteroides_stools <- genus_over_time("Bacteroides", genus_stools_t)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_bacteroides_stools.png", bacteroides_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
bacteroides_tissues <- genus_over_time("Bacteroides", genus_tissues_t)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_bacteroides_tissues.png", bacteroides_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Clostridiales over time
clostridiales_stools <- genus_over_time("Clostridiales Unclassified", genus_stools_t)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_clostridiales_stools.png", clostridiales_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
clostridiales_tissues <- genus_over_time("Clostridiales Unclassified", genus_tissues_t)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_clostridiales_tissues.png", clostridiales_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Ruminococcaceae over time
ruminococcaceae_stools <- genus_over_time("Ruminococcaceae Unclassified", genus_stools_t)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_ruminococcaceae_stools.png", ruminococcaceae_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
ruminococcaceae_tissues <- genus_over_time("Ruminococcaceae Unclassified", genus_tissues_t)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_ruminococcaceae_tissues.png", ruminococcaceae_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Peptostreptococcaceae over time
peptostreptococcaceae_stools <- genus_over_time("Peptostreptococcaceae Unclassified", genus_stools_t)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_peptostreptococcaceae_stools.png", peptostreptococcaceae_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
peptostreptococcaceae_tissues <- genus_over_time("Peptostreptococcaceae Unclassified", genus_tissues_t)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_peptostreptococcaceae_tissues.png", peptostreptococcaceae_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Acetatifactor over time
acetatifactor_stools <- genus_over_time("Acetatifactor", genus_stools_t)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_acetatifactor_stools.png", acetatifactor_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
acetatifactor_tissues <- genus_over_time("Acetatifactor", genus_tissues_t)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_acetatifactor_tissues.png", acetatifactor_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Akkermansia over time
akkermansia_stools <- genus_over_time("Akkermansia", genus_stools_t)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_akkermansia_stools.png", akkermansia_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
akkermansia_tissues <- genus_over_time("Akkermansia", genus_tissues_t)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_genus_akkermansia_tissues.png", akkermansia_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Analysis of mock-infected mice compared to their corresponding group that was challenged with C. difficile----
#Restrict mock analysis to C, WM, CN, and WMN mice. 

#Custom scale for mock plots----
color_scheme_m <- c("#238b45", "#88419d", "#238b45", "#88419d")
color_groups_m <- c("CN", "WMN", "C", "WM")
color_labels_m <- c("Clind. without infection", "5-day PEG without infection", "Clind.", "5-day PEG")
linetype_scheme_m <- c("longdash", "solid")
linetype_infected <- c("no", "yes")

#Also create dataframes of diversity data that includes mock challenged mice (WMN and C), separated into stool and tissue samples
diversity_mock_stools <- subset_stool(add_mocks(diversity_subset, diversity_data)) %>% 
  filter(group %in% c("C", "CN", "WM", "WMN"))
diversity_mock_tissues <- subset_tissue(add_mocks(diversity_subset, diversity_data)) %>% 
  filter(group %in% c("C", "CN", "WM", "WMN"))

num_mock_stool <- count_subset(diversity_mock_stools)#Number of stool samples per group per day + mock challenged mice
stool_mock_test_days <- c(-5, -1, 0, 4, 6, 30)

num_mock_tissue <- count_subset(diversity_mock_tissues) #Number of stool samples per group per day + mock challenged mice
tissue_mock_test_days <- c(4, 6, 30)

#Function to perform Kruskal-Wallis test for differences in Shannon diversity index across infected and mock groups on a particular day with Benjamini Hochberg correction
#Function will then perform pairwise comparisons for any days where there was a difference between groups
#Arguments: 
#diversity_subset <- subset (stools or tissue samples) of diversity_data to perform statistical test on
#timepoint = timepoints to assess differences between groups specific to the subset (stool or tissue)
#subset_name = label to append to results filename to indicate subset analyzed. Ex. stool, tissues, stool_mock, tissues_mock
stats_shannon_mock <- function(diversity_subset, timepoints, subset_name){
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
    select(-data) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) 
  #Pull significant days
  diversity_sig_days <- pull_sig_days(diversity_stats_adjust)
  diversity_stats_pairwise <- diversity_stats %>% 
    filter(day %in% diversity_sig_days) %>% #only perform pairwise tests for days that were significant 
    group_by(day) %>% 
    mutate(model=map(data, ~pairwise.wilcox.test(x=.x$shannon, g=as.factor(.x$group), p.adjust.method="BH") %>% 
                       tidy() %>% 
                       mutate(compare=paste(group1, group2, sep="-")) %>% 
                       select(-group1, -group2) %>% 
                       pivot_wider(names_from=compare, values_from=p.value)
    )
    ) %>% 
    unnest(model) %>% 
    select(-data, -parameter, -statistic, -p.value, -method, -WM, -C, -CN, -WMN) %>% 
    #Combine with diversity_stats_adjust so that adjusted p-values are on the same table
    inner_join(diversity_stats_adjust, by = c("day")) %>% 
    select(-p.value, -parameter, -statistic) %>% 
    write_tsv(path = paste0("data/process/5_days_PEG_shannon_stats_", subset_name, "_subset.tsv"))
}

#Test with shannon for mock stool subset
mock_kw_shannon_stools <- stats_shannon_mock(diversity_mock_stools, stool_mock_test_days, "mock_stools")

#Function to perform Kruskal-Wallis test for differences in Shannon diversity index across infected and mock groups on a particular day with Benjamini Hochberg correction
#Function will then perform pairwise comparisons for any days where there was a difference between groups
#Arguments: 
#diversity_subset <- subset (stools or tissue samples) of diversity_data to perform statistical test on
#timepoint = timepoints to assess differences between groups specific to the subset (stool or tissue)
#subset_name = label to append to results filename to indicate subset analyzed. Ex. stool, tissues, stool_mock, tissues_mock
stats_shannon_mock_t <- function(diversity_subset, timepoints, subset_name){
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
    select(-data) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) 
  #Pull significant days
  diversity_sig_days <- pull_sig_days(diversity_stats_adjust)
  diversity_stats_pairwise <- diversity_stats %>% 
    filter(day %in% diversity_sig_days) %>% #only perform pairwise tests for days that were significant 
    group_by(day) %>% 
    mutate(model=map(data, ~pairwise.wilcox.test(x=.x$shannon, g=as.factor(.x$group), p.adjust.method="BH") %>% 
                       tidy() %>% 
                       mutate(compare=paste(group1, group2, sep="-")) %>% 
                       select(-group1, -group2) %>% 
                       pivot_wider(names_from=compare, values_from=p.value)
    )
    ) %>% 
    unnest(model) %>% 
    select(-data, -parameter, -statistic, -p.value, -method, -WM, -C, -WMN) %>% 
    #Combine with diversity_stats_adjust so that adjusted p-values are on the same table
    inner_join(diversity_stats_adjust, by = c("day")) %>% 
    select(-p.value, -parameter, -statistic) %>% 
    write_tsv(path = paste0("data/process/5_days_PEG_shannon_stats_", subset_name, "_subset.tsv"))
}
#Test with shannon for mock tissue subset
mock_kw_shannon_tissues <- stats_shannon_mock_t(diversity_mock_tissues, tissue_test_days, "mock_tissues")

#Plot of shannon diversity over time
shannon_mock_stools_plot <- diversity_mock_stools %>% 
  mutate(day = as.integer(day)) %>% 
  group_by(group, day) %>%
  mutate(median_shannon = median(shannon)) %>%
  ggplot(x = day, y = shannon, colour = group)+
  geom_line(mapping = aes(x = day, y = median_shannon, group = group, color = group, linetype = infected), alpha = 0.6, size = 1) +
  scale_colour_manual(name=NULL,
                      values=color_scheme_m,
                      breaks=color_groups_m,
                      labels=color_labels_m) +
  scale_linetype_manual(name="Infected", #Dashed lines are mock challenged mice
                        values=linetype_scheme_m,
                        breaks=linetype_infected,
                        labels=linetype_infected) +
  scale_y_continuous(limits = c(0,4.1), expand = c(0, 0))+ #expand argument gets rid of the extra space around the scale
  theme_classic()+
  labs(title=NULL,
       x="Days Post-Infection",
       y="Shannon Diversity Index")+
  theme(legend.position = "none", #Remove legend
        text = element_text(size = 16), # Change font size for entire plot
        axis.ticks.x = element_blank())+
  scale_x_continuous(breaks = c(-5, -1, 0, 4, 6, 30),
                     limits = c(-5.5,31))
save_plot(filename = "results/figures/5_days_PEG_shannon_stools_mock.png", shannon_mock_stools_plot, base_height = 4, base_width = 9, base_aspect_ratio = 2)
#Plot of shannon diversity in mock tissues over time
shannon_mock_tissues_plot <- diversity_mock_tissues %>% 
  mutate(day = as.integer(day)) %>% 
  group_by(group, day) %>%
  mutate(median_shannon = median(shannon)) %>%
  ggplot(x = day, y = shannon, colour = group)+
  geom_line(mapping = aes(x = day, y = median_shannon, group = group, color = group, linetype = infected), alpha = 0.6, size = 1) +
  scale_colour_manual(name=NULL,
                      values=color_scheme_m,
                      breaks=color_groups_m,
                      labels=color_labels_m) +
  scale_linetype_manual(name="Infected", #Dashed lines are mock challenged mice
                        values=linetype_scheme_m,
                        breaks=linetype_infected,
                        labels=linetype_infected) +
  scale_y_continuous(limits = c(0,4.1), expand = c(0, 0))+ #expand argument gets rid of the extra space around the scale
  theme_classic()+
  labs(title=NULL,
       x="Days Post-Infection",
       y="Shannon Diversity Index")+
  theme(legend.position = "none", #Remove legend
        text = element_text(size = 16), # Change font size for entire plot
        axis.ticks.x = element_blank())+
  scale_x_continuous(breaks = c(0, 4, 6, 30),
                     limits = c(-0.5,31))
save_plot(filename = "results/figures/5_days_PEG_shannon_tissues_mock.png", shannon_mock_tissues_plot, base_height = 4, base_width = 9, base_aspect_ratio = 2)

#Function to perform Kruskal-Wallis test for differences in richness diversity index across infected and mock groups on a particular day with Benjamini Hochberg correction
#Function will then perform pairwise comparisons for any days where there was a difference between groups
#Arguments: 
#diversity_subset <- subset (stools or tissue samples) of diversity_data to perform statistical test on
#timepoint = timepoints to assess differences between groups specific to the subset (stool or tissue)
#subset_name = label to append to results filename to indicate subset analyzed. Ex. stool, tissues, stool_mock, tissues_mock
stats_richness_mock <- function(diversity_subset, timepoints, subset_name){
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
    select(-data) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) 
  #Pull significant days
  diversity_sig_days <- pull_sig_days(diversity_stats_adjust)
  diversity_stats_pairwise <- diversity_stats %>% 
    filter(day %in% diversity_sig_days) %>% #only perform pairwise tests for days that were significant 
    group_by(day) %>% 
    mutate(model=map(data, ~pairwise.wilcox.test(x=.x$sobs, g=as.factor(.x$group), p.adjust.method="BH") %>% 
                       tidy() %>% 
                       mutate(compare=paste(group1, group2, sep="-")) %>% 
                       select(-group1, -group2) %>% 
                       pivot_wider(names_from=compare, values_from=p.value)
    )
    ) %>% 
    unnest(model) %>% 
    select(-data, -parameter, -statistic, -p.value, -method, -WM, -C, -CN, -WMN) %>% 
    #Combine with diversity_stats_adjust so that adjusted p-values are on the same table
    inner_join(diversity_stats_adjust, by = c("day")) %>% 
    select(-p.value, -parameter, -statistic) %>% 
    write_tsv(path = paste0("data/process/5_days_PEG_richness_stats_", subset_name, "_subset.tsv"))
}

#Test with richness for mock stool subset
mock_kw_richness_stools <- stats_richness_mock(diversity_mock_stools, stool_mock_test_days, "mock_stools")

#Function to perform Kruskal-Wallis test for differences in richness diversity index across infected and mock groups on a particular day with Benjamini Hochberg correction
#Function will then perform pairwise comparisons for any days where there was a difference between groups
#Arguments: 
#diversity_subset <- subset (stools or tissue samples) of diversity_data to perform statistical test on
#timepoint = timepoints to assess differences between groups specific to the subset (stool or tissue)
#subset_name = label to append to results filename to indicate subset analyzed. Ex. stool, tissues, stool_mock, tissues_mock
stats_richness_mock_t <- function(diversity_subset, timepoints, subset_name){
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
    select(-data) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) 
  #Pull significant days
  diversity_sig_days <- pull_sig_days(diversity_stats_adjust)
  diversity_stats_pairwise <- diversity_stats %>% 
    filter(day %in% diversity_sig_days) %>% #only perform pairwise tests for days that were significant 
    group_by(day) %>% 
    mutate(model=map(data, ~pairwise.wilcox.test(x=.x$sobs, g=as.factor(.x$group), p.adjust.method="BH") %>% 
                       tidy() %>% 
                       mutate(compare=paste(group1, group2, sep="-")) %>% 
                       select(-group1, -group2) %>% 
                       pivot_wider(names_from=compare, values_from=p.value)
    )
    ) %>% 
    unnest(model) %>% 
    select(-data, -parameter, -statistic, -p.value, -method, -WM, -C, -WMN) %>% 
    #Combine with diversity_stats_adjust so that adjusted p-values are on the same table
    inner_join(diversity_stats_adjust, by = c("day")) %>% 
    select(-p.value, -parameter, -statistic) %>% 
    write_tsv(path = paste0("data/process/5_days_PEG_richness_stats_", subset_name, "_subset.tsv"))
}

#Test with richness for mock tissue subset
mock_kw_richness_tissues <- stats_richness_mock_t(diversity_mock_tissues, tissue_test_days, "mock_tissues")

#Plot of richness over time
richness_mock_stools_plot <- diversity_mock_stools %>% 
  mutate(day = as.integer(day)) %>% 
  group_by(group, day) %>%
  mutate(median_sobs = median(sobs)) %>%
  ggplot(x = day, y = sobs, colour = group)+
  geom_line(mapping = aes(x = day, y = median_sobs, group = group, color = group, linetype = infected), alpha = 0.6, size = 1) +
  scale_colour_manual(name=NULL,
                      values=color_scheme_m,
                      breaks=color_groups_m,
                      labels=color_labels_m) +
  scale_linetype_manual(name="Infected", #Dashed lines are mock challenged mice
                        values=linetype_scheme_m,
                        breaks=linetype_infected,
                        labels=linetype_infected) +
  scale_y_continuous(limits = c(0,160), expand = c(0, 0))+ #expand argument gets rid of the extra space around the scale
  theme_classic()+
  labs(title=NULL,
       x="Days Post-Infection",
       y="Number of Observed OTUs")+
  theme(legend.position = "none", #Remove legend
        text = element_text(size = 16), # Change font size for entire plot
        axis.ticks.x = element_blank())+
  scale_x_continuous(breaks = c(-5, -1, 0, 4, 6, 30),
                     limits = c(-5.5,31))
save_plot(filename = "results/figures/5_days_PEG_richness_stools_mock.png", richness_mock_stools_plot, base_height = 4, base_width = 9, base_aspect_ratio = 2)

richness_mock_tissues_plot <- diversity_mock_tissues %>% 
  mutate(day = as.integer(day)) %>% 
  group_by(group, day) %>%
  mutate(median_sobs = median(sobs)) %>%
  ggplot(x = day, y = sobs, colour = group)+
  geom_line(mapping = aes(x = day, y = median_sobs, group = group, color = group, linetype = infected), alpha = 0.6, size = 1) +
  scale_colour_manual(name=NULL,
                      values=color_scheme_m,
                      breaks=color_groups_m,
                      labels=color_labels_m) +
  scale_linetype_manual(name="Infected", #Dashed lines are mock challenged mice
                        values=linetype_scheme_m,
                        breaks=linetype_infected,
                        labels=linetype_infected) +
  scale_y_continuous(limits = c(0,160), expand = c(0, 0))+ #expand argument gets rid of the extra space around the scale
  theme_classic()+
  labs(title=NULL,
       x="Days Post-Infection",
       y="Number of Observed OTUs")+
  theme(legend.position = "none", #Remove legend
        text = element_text(size = 16), # Change font size for entire plot
        axis.ticks.x = element_blank())+
  scale_x_continuous(breaks = c(0, 4, 6, 30),
                     limits = c(-0.5,31))
save_plot(filename = "results/figures/5_days_PEG_richness_tissues_mock.png", richness_mock_tissues_plot, base_height = 4, base_width = 9, base_aspect_ratio = 2)

#PCoA analysis of mock challenged mice----
#Define shape scheme based on Infected status
shape_scheme <- c(1, 19)
shape_infected <- c("no", "yes")

#Function to plot PCoA data for mock-challenged mice comparison
plot_mock_pcoa <- function(df){
  ggplot(df, aes(x=axis1, y=axis2, color = group, alpha = day, shape = infected)) +
    geom_point(size=2) +
    scale_colour_manual(name=NULL,
                        values=color_scheme_m,
                        breaks=color_groups_m,
                        labels=color_labels_m)+ 
    scale_shape_manual(name="Infected",
                       values=shape_scheme,
                       breaks=shape_infected,
                       labels=shape_infected) +
    coord_fixed() +
    labs(x="PCoA 1",
         y="PCoA 2",
         alpha= "Day") +
    theme_classic()+
    theme(legend.position = "bottom",
          text = element_text(size = 16))
}

pcoa_mock_stool <- read_tsv("data/process/5_day_PEG/stools_mock/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  mutate(day = as.integer(day)) %>% #Transform day into continuous
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

pcoa_axes_5_day_PEG_stool <- read_tsv("data/process/5_day_PEG/stools_mock/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

plot_pcoa_mock_stool <- plot_mock_pcoa(pcoa_mock_stool)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  theme(legend.position = "none")
save_plot(filename = paste0("results/figures/5_days_PEG_stool_PCoA_mock.png"), plot_pcoa_mock_stool, base_height = 5, base_width = 5)

#Create stand alone legend for mock PCoAs
mock_group_legend <- pcoa_mock_stool  %>%
  filter(group %in% c("WM", "C")) %>% #Just need these 2 groups for color legend
  ggplot(aes(x = axis1, y = axis2, color=group, alpha = day))+
  scale_colour_manual(name=NULL,
                      values=color_scheme_m,
                      breaks=color_groups_m,
                      labels=color_labels_m)+ 
  geom_point()+ theme_classic()+ theme(legend.position = "bottom")
mock_group_legend <- get_legend(mock_group_legend)
mock_shape_legend <- pcoa_mock_stool  %>%
  ggplot(aes(x = axis1, y = axis2, shape = infected,))+
  scale_shape_manual(name="Infected",
                     values=shape_scheme,
                     breaks=shape_infected,
                     labels=shape_infected) +
  geom_point()+ theme_classic()+ theme(legend.position = "bottom")
mock_shape_legend <- get_legend(mock_shape_legend)
mock_legend <- plot_grid(mock_group_legend, mock_shape_legend, nrow = 2)
save_plot("results/figures/5_days_PEG_pcoa_mock_legend.png", mock_legend, base_height = 1, base_width = 5)

#Create PCoA of mock tisuse samples
pcoa_mock_tissue <- read_tsv("data/process/5_day_PEG/tissues_mock/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  mutate(day = as.integer(day)) %>% #Transform day into continuous
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

pcoa_axes_5_day_PEG_tissue <- read_tsv("data/process/5_day_PEG/tissues_mock/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

plot_pcoa_mock_tissue <- plot_mock_pcoa(pcoa_mock_tissue)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  theme(legend.position = "none")
save_plot(filename = paste0("results/figures/5_days_PEG_tissue_PCoA_mock.png"), plot_pcoa_mock_tissue, base_height = 5, base_width = 5)

#Create stand alone alpha legend for mock PCoA of tissue samples
mock_alpha_legend <- pcoa_mock_tissue  %>%
  ggplot(aes(x = axis1, y = axis2, alpha = day))+
  geom_point()+ theme_classic()+ theme(legend.position = "bottom")
mock_alpha_legend <- get_legend(mock_alpha_legend)
save_plot("results/figures/5_days_PEG_pcoa_mock_alpha_legend_tissues.png", mock_alpha_legend, base_height = 0.5, base_width = 2.5)

#Line Plots of specific bacterial genera in mock-infected mice----

#Function to create faceted line plots of relative abundances of genera of interest over time for mock challenged mice comparisons
#sample_df = subset dataframe of samples to be plotted
#genera = list of genera to plot
#timepoints = days of the experiment to plot
line_plot_mock_genus <- function(sample_df, genera, timepoints){
  sample_df %>% 
    filter(group %in% c("C", "WM", "CN","WMN")) %>%  ##Select only mock challenged & corresponding C. diff challenged groups
    filter(genus %in% genera) %>% #Select only genera of interest
    mutate(genus = fct_relevel(genus, genera)) %>% #Reorder genera to match order of genera of interest
    mutate(group = fct_relevel(group, rev(color_groups))) %>% #Specify the order of the groups
    filter(day %in% timepoints) %>% #Select only timepoints of interest
    mutate(day = as.integer(day)) %>% 
    group_by(group, genus, day, infected) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_line(aes(x = day, y=median, color=group, linetype =infected), alpha = 0.6, size = 1)+
    scale_colour_manual(name=NULL,
                        values=color_scheme_m,
                        breaks=color_groups_m,
                        labels=color_labels_m) +
    scale_linetype_manual(name=NULL, #Dashed lines are mock challenged mice
                          values=linetype_scheme_m,
                          breaks=linetype_infected,
                          labels=linetype_infected) +
    scale_y_continuous(trans = "log10", limits = c(1/10900, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    geom_hline(yintercept=1/1000, color="gray")+ #Represents limit of detection
    labs(title=NULL,
         x="Days Post-Infection",
         y="Relative abundance (%)")+
    facet_wrap(~genus, nrow = 2, labeller = label_wrap_gen(width = 12))+
    theme_classic()+
    theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
          strip.text = element_text(face = "italic"),
          plot.title = element_markdown(hjust = 0.5), #Have only the genera names show up as italics
          text = element_text(size = 16),
          legend.position = "None")
}

#Lineplots of the top 6 significant genera in mock stool samples----
lp_stool_days <- c("-5","-1", "0", "4", "6", "30")
facet_labels <- top_10_sig_genus[1:6] #Pick just the top 6
names(facet_labels) <- top_10_sig_genus[1:6] #Pick just the top 6
lp_stool_mock <- line_plot_mock_genus(genus_mock_stools, 
                                 top_10_sig_genus[1:6], lp_stool_days)+
  scale_x_continuous(limits = c(-5.5,30.5), breaks = c(-5, 0, 4, 6, 30), labels = c(-5, 0, 4, 6, 30))#Change scale
save_plot(filename = "results/figures/5_days_PEG_genus_lineplot_mock_stools.png", lp_stool_mock, base_height = 5, base_width =10)

#Lineplots of the top 6 significant genera in mock tissue samples----
lp_tissue_days <- diversity_tissues %>% distinct(day) %>% 
  filter(day %in% c("4", "6", "30")) %>% #Focus on day -1 through 10 timepoints
  pull(day)
facet_labels <- top_sig_genus_tissues[1:6] #Pick just the top 6
names(facet_labels) <- top_sig_genus_tissues[1:6] #Pick just the top 6
lp_tissue_mock <- line_plot_mock_genus(genus_mock_tissues, 
                                    top_sig_genus_tissues[1:6], lp_tissue_days)+
  scale_x_continuous(limits = c(3.5, 30.5), breaks = c(4, 6, 30), labels = c(4, 6, 30))#Change scale since we have less timepoints for tissues
save_plot(filename = "results/figures/5_days_PEG_genus_lineplot_mock_tissues.png", lp_tissue_mock, base_height = 5, base_width = 8)

