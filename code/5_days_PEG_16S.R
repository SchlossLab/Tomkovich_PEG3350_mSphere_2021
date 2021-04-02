source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match 5_days_PEG plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMR")
color_labels <- c( "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 + 10-day recovery")
#Need to create an additional color scheme with 6 colors (or consider keeping colors and doing open/closed for mock challenged mice)
#See 5_days_PEG_histology_scores.R for how mock challenged mice were presented
#Define shape scheme based on Infected status----
shape_scheme <- c(1, 19)
shape_infected <- c("no", "yes")

#metadata <- metadata %>%
#  mutate(day = as.integer(day))  #Day variable (transformed to integer to get rid of decimals on PCoA animation
#agg_otu_data <- agg_otu_data %>% 
#  mutate(day = as.integer(day))  #Day variable (transformed to integer to get rid of decimals on PCoA animation

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
#Also create dataframes of diversity data that includes mock challenged mice (WMN and C), separated into stool and tissue samples
diversity_mock_stools <- subset_stool(add_mocks(diversity_subset, diversity_data))
diversity_mock_tissues <- subset_tissue(add_mocks(diversity_subset, diversity_data))

#Figure out how many samples we have per group per day for each subset
num_stool <- count_subset(diversity_stools) #Number of stool samples per group per day
num_tissue <- count_subset(diversity_tissues) #Number of tissue samples per group per day
num_mock_stool <- count_subset(diversity_mock_stools)#Number of stool samples per group per day + mock challenged mice
num_mock_tissue <- count_subset(diversity_mock_tissues) #Number of stool samples per group per day + mock challenged mice

#Experimental days to analyze with the Kruskal-Wallis test (timepoints with 16S data for at least 3 groups)
#Baseline (before treatment) for WMR is day -15. For C, WM, and WMC baseline is day -5
stool_test_days <- c(-5, -1, 0, 1, 2, 3, 4, 5, 6, 10, 30)
#stool_mock_test_days #Don't test, included in code/all_subsets_16S.R
tissue_test_days <- c(6, 30) #Only 2 days with samples from at least 3 groups
#tissue_mock_test_days #Don't test, included in code/all_subsets_16S.R

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
    select(day, statistic, p.value, parameter, method, C, WM, WMC, WMR) %>% 
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
    select(day, statistic, p.value, parameter, method, C, WM, WMC, WMR) %>% 
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
save_plot(filename = "results/figures/5_days_PEG_shannon_stools.png", shannon_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

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

#Read in pcoa loadings and axes for 5_day_PEG PCoA subset
#Pull 5_Day_PEG subset of PCoA data
pcoa_5_day_PEG <- read_tsv("data/process/5_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>% 
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and use left_join to keep all samples in pcoa data frame
  filter(!is.na(axis1)) %>%  #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff
  filter(!group %in% c("WMN", "CN")) #Remove the mock challenged mice
  
#Pull axes from loadings file
pcoa_axes_5_day_PEG <- read_tsv("data/process/5_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot <- plot_pcoa(pcoa_5_day_PEG)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))

save_plot(filename = paste0("results/figures/5_Day_PEG_PCoA.png"), pcoa_subset_plot, base_height = 5, base_width = 5)

#Create stand alone legend
group_legend <- pcoa_5_day_PEG  %>%
  ggplot(aes(x = axis1, y = axis2, color = group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  geom_point()+ theme_classic()
group_legend <- get_legend(group_legend)
save_plot("results/figures/5_days_PEG_pcoa_legend.png", group_legend, base_height = 1, base_width = 2.3)

#Pull 5_Day_PEG subset of PCoA 
#Stool Subset
pcoa_5_day_PEG_stool <- read_tsv("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

pcoa_5_day_PEG_stool<- subset(pcoa_5_day_PEG_stool, !group %in% c("WMN", "CN"))

pcoa_axes_5_day_PEG_stool <- read_tsv("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot_stool <- plot_pcoa(pcoa_5_day_PEG_stool)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))

save_plot(filename = paste0("results/figures/5_Day_PEG_stool_PCoA.png"), pcoa_subset_plot_stool, base_height = 5, base_width = 5)

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

#Tissue subset
pcoa_5_day_PEG_tissues <- read_tsv("data/process/5_day_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and use left_join to keep all samples in pcoa data frame
  filter(!is.na(axis1)) %>%  #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff
  filter(!group %in% c("WMN", "CN")) #Remove the mock challenged mice

#Pull axes from loadings file
pcoa_axes_5_day_PEG_tissues <- read_tsv("data/process/5_day_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG_tissues %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG_tissues %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot_tissue <- plot_pcoa(pcoa_5_day_PEG_tissues)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  scale_alpha_continuous(range = c(.3, 1),
                         breaks= c(4, 6, 20, 30),
                         labels=c(4, 6, 20, 30))

save_plot(filename = paste0("results/figures/5_Day_PEG_tissues_PCoA.png"), pcoa_subset_plot_tissue, base_height = 5, base_width = 5)

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

#Function to test for differences across groups at the OTU level for specific timepoints
#Function to test at the otu level:
#Arguments:
# timepoint = day of the experiment
#sample_df = subset dataframe of just stool or tissue samples
#sample_type = "stool" or "tissue" to be included in filename
kruskal_wallis_otu <- function(timepoint, sample_df, sample_type){
  otu_stats <- sample_df %>%
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
    select(-data) %>% #Keep everything but the data column
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
    arrange(p.value.adj) %>%
    write_tsv(path = paste0("data/process/5_days_PEG_otu_stats_day_", timepoint, "_", sample_type, ".tsv"))
}

#Create empty data frame to combine stat dataframes for all days that were tested
kw_otu_stools <- data.frame(otu=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                                     WM =double(),C =double(),WMR =double(),WMC=double(),
                            p.value.adj=double(),day=double())

# Perform kruskal wallis tests at the otu level for the stool samples----
for (d in stool_test_days){
  kruskal_wallis_otu(d, otu_stools, "stools")
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/5_days_PEG_otu_stats_day_", d, "_stools.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_otu_stools_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, otu))
  kw_otu_stools <- add_row(kw_otu_stools, stats)  #combine all the dataframes together
}

#Create empty data frame to combine stat dataframes for all days that were tested
kw_otu_tissues <- data.frame(otu=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                            WM =double(),C =double(),WMR =double(),WMC=double(),
                            p.value.adj=double(),day=double())
# Perform kruskal wallis tests at the otu level for the tissue samples----
for (d in tissue_test_days){
  kruskal_wallis_otu(d, otu_tissues, "tissues")
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/5_days_PEG_otu_stats_day_", d, "_tissues.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_otu_tissues_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, otu))
  kw_otu_tissues <- add_row(kw_otu_tissues, stats)  #combine all the dataframes together
}

#OTUs that varied across treatment groups and were shared across days 
#Stools
view(sig_otu_stools_day1)
view(sig_otu_stools_day10)
view(sig_otu_stools_day4)
shared_sig_stools_otus_D1toD6 <- intersect_all(sig_otu_stools_day1, sig_otu_stools_day2,                                              sig_otu_stools_day3, sig_otu_stools_day4, 
                                             sig_otu_stools_day5, sig_otu_tissues_day6) #fill in different days to compare
view(shared_sig_stools_otus_D1toD6)
print(shared_sig_stools_otus_D1toD6)

# "Bacteroides (OTU 1)"            "Peptostreptococcaceae (OTU 12)"
# "Enterobacteriaceae (OTU 2)"     "Lachnospiraceae (OTU 20)"      
# "Lachnospiraceae (OTU 4)"        "Lachnospiraceae (OTU 32)"


#Tissues
shared_sig_tissues_otus <- intersect_all(sig_otu_tissues_day30, sig_otu_tissues_day6) #fill in different days to compare
view(shared_sig_tissues_otus)

#Plots of the relative abundances of OTUs that significantly varied across sources of mice from day -1 to day 1----
otus_d1 <- plot_otus_dx(otu_stools, sig_otu_stools_day1[1:20], 1)+
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  ggtitle("Day 1 post-infection Stools")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_otus_stools_d1_top20.png", otus_d1, base_height = 7, base_width = 8)

otus_d4 <- plot_otus_dx(otu_stools, sig_otu_stools_day4 [1:20], 4)+
  ggtitle("Day 4 post-infection Stools")+ #Title plot
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_otus_stools_d4_top20.png", otus_d4, base_height = 7, base_width = 8)

otus_d10 <- plot_otus_dx(otu_stools, sig_otu_stools_day10[1:20], 10)+
  ggtitle("Day 10 post-infection Stools")+ #Title plot
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_otus_stools_d10.png", otus_d10, base_height = 7, base_width = 8)

otus_d6 <- plot_otus_dx(otu_stools, `sig_otu_stools_day6` [1:20], 6)+
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  ggtitle("Day 6 post-infection Stools")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_otus_stools_d6_top20.png", otus_d6, base_height = 7, base_width = 8)

otus_tissues_d6 <- plot_otus_dx(otu_tissues, sig_otu_tissues_day6[1:20], 6)+
  ggtitle("Day 6 post-infection Tissue")+ #Title plot
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_otus_tissues_d6.png", otus_tissues_d6, base_height = 7, base_width = 8)

otus_tissues_d30 <- plot_otus_dx(otu_tissues, sig_otu_tissues_day30[1:20], 30)+
  ggtitle("Day 30 post-infection Tissue")+ #Title plot
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_otus_tissues_d30.png", otus_tissues_d30, base_height = 7, base_width = 8)

#Heatmaps of significant OTUs over time facet by group----
#Create list of otus to plot
all_sig_otus <- c(`sig_otu_stools_day-5`, `sig_otu_stools_day-1`, sig_otu_stools_day0, sig_otu_stools_day1, 
                  `sig_otu_stools_day2`, `sig_otu_stools_day3`, sig_otu_stools_day4, sig_otu_stools_day5,
                  `sig_otu_stools_day6`, `sig_otu_stools_day10`, sig_otu_stools_day30)
#370 total OTUs
unique_sig_otus <-unique(all_sig_otus)
#130 unique significant OTUs
#Rank the 130 signifiant OTUs in order of agg_rel_abund, select the top 20
hm_sig_otus_abund <-  otu_stools %>% 
  filter(otu %in% unique_sig_otus) %>%
  group_by(otu) %>% 
  summarize(median=(median(agg_rel_abund + 1/2000))) %>% #Median relative abundance for all samples for a particular OTU
  arrange(desc(median)) %>% #Arrange largest to smallest
  slice_head(n = 20) %>% 
  pull(otu)

#Rank OTUs by adjusted p-value
hm_sig_otus_p_adj <- kw_otu_stools %>% 
  filter(p.value.adj < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(otu) %>% 
  slice_head(n = 25) %>% 
  pull(otu)

hm_stool_days <- diversity_stools %>% distinct(day) %>% pull(day)
facet_labels <- color_labels #Create descriptive labels for facets
names(facet_labels) <- c("C", "WM", "WMC", "WMR") #values that correspond to group, which is the variable we're faceting by
hm_stool <- hm_plot_otus(otu_stools, hm_sig_otus_p_adj, hm_stool_days)+
  scale_x_discrete(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30), labels = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_otus_heatmap_stools.png", hm_stool, base_height = 14, base_width = 15)

#Create heatmap of significant OTUs for tissue samples----
#Rank OTUs by adjusted p-value
hm_sig_otus_p_adj_tissues <- kw_otu_tissues %>% 
  filter(p.value.adj < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(otu) %>% 
  slice_head(n = 25) %>% 
  pull(otu)
hm_tissue_days <- diversity_tissues %>% distinct(day) %>% pull(day)
hm_tissues <- hm_plot_otus(otu_tissues, hm_sig_otus_p_adj_tissues, hm_tissue_days)+
  scale_x_discrete(breaks = c(4, 6, 20, 30), labels = c(4, 6, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_otus_heatmap_tissues.png", hm_tissues, base_height = 14, base_width = 15)
#Examine OTUs that were significant in stool samples
hm_tissues_stool_otus <- hm_plot_otus(otu_tissues, hm_sig_otus_p_adj, hm_tissue_days)+
  scale_x_discrete(breaks = c(4, 6, 20, 30), labels = c(4, 6, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_otus_heatmap_tissues_stool_otus.png", hm_tissues_stool_otus, base_height = 14, base_width = 15)
#Pull list of otus that overlap between stool & tissue heatmaps
hm_overlap <- intersect_all(hm_sig_otus_p_adj, hm_sig_otus_p_adj_tissues)
#11 OTUs overlap:  
#"Lachnospiraceae (OTU 33)"       "Blautia (OTU 19)"              
#"Ruminococcaceae (OTU 50)"       "Ruminococcaceae (OTU 54)"      
#"Ruminococcaceae (OTU 92)"       "Oscillibacter (OTU 45)"        
# "Lachnospiraceae (OTU 30)"       "Lachnospiraceae (OTU 31)"      
# "Lachnospiraceae (OTU 4)"        "Peptostreptococcaceae (OTU 12)"
# "Enterobacteriaceae (OTU 2)" 

#Create heatmaps of mock challenged mice----
#Mock stool samples
otu_mock_only_stools  <-  otu_mock_stools %>% 
  filter(group %in% c("WMN", "CN"))
hm_mock_stool_days <- otu_mock_only_stools %>% distinct(day) %>% pull(day)
facet_labels <- c("5-day PEG 3350 without infection", "Clind. without infection") #Create descriptive labels for facets
names(facet_labels) <- c("WMN", "CN") #values that correspond to group, which is the variable we're faceting by
hm_mock_stool <- hm_plot_otus(otu_mock_only_stools, hm_sig_otus_p_adj, hm_mock_stool_days)+
  scale_x_discrete(breaks = c(-5, -1, 0, 4, 6, 30), labels = c(-5, -1, 0, 4, 6, 30)) 
save_plot(filename = "results/figures/5_days_PEG_otus_heatmap_stools_mock.png", hm_mock_stool, base_height = 14, base_width = 15)

#Mock tissue samples
otu_mock_only_tissues <- otu_mock_tissues %>% 
  filter(group %in% c("WMN", "CN")) 
hm_mock_tissue_days <- otu_mock_only_tissues %>% distinct(day) %>% pull(day)
hm_mock_tissues <- hm_plot_otus(otu_mock_only_tissues, hm_sig_otus_p_adj_tissues, hm_mock_tissue_days)+
  scale_x_discrete(breaks = c(0, 4, 6, 30), labels = c(0, 4, 6, 30)) 
save_plot(filename = "results/figures/5_days_PEG_otus_heatmap_tissues_mock.png", hm_mock_tissues, base_height = 14, base_width = 15)



#Examine C. difficile otu over time----
peptostrep_stools <- otu_over_time("Peptostreptococcaceae (OTU 12)", otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_peptostreptococcaceae_stools.png", peptostrep_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
peptostrep_tissues <- otu_over_time("Peptostreptococcaceae (OTU 12)", otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_peptostreptococcaceae_tissues.png", peptostrep_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Bacteroides otu over time----
bacteroides_stools <- otu_over_time("Bacteroides (OTU 1)", otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_bacteroides_stools.png", bacteroides_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
bacteroides_tissues <- otu_over_time("Bacteroides (OTU 1)", otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_bacteroides_tissues.png", bacteroides_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Enterobacteriaceae (OTU 2) over time----
entero2_stools <- otu_over_time("Enterobacteriaceae (OTU 2)", otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_entero2_stools.png", entero2_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
entero2_tissues <- otu_over_time("Enterobacteriaceae (OTU 2)", otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_entero2_tissues.png", entero2_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine Porphyromonadaceae otu over time----
porph_stools <- otu_over_time("Porphyromonadaceae (OTU 8)", otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_porph_stools.png", porph_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
porph_tissues <- otu_over_time("Porphyromonadaceae (OTU 8)", otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_porph_tissues.png", porph_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine "Lachnospiraceae (OTU 20)"  over time----
lach20_stools <- otu_over_time("Lachnospiraceae (OTU 20)" , otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_lach20_stools.png", lach20_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
lach20_tissues <- otu_over_time("Lachnospiraceae (OTU 20)" , otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_lach20_tissues.png", lach20_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Examine "Lachnospiraceae (OTU 4)"  over time----
lach4_stools <- otu_over_time("Lachnospiraceae (OTU 4)" , otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_lach4_stools.png", lach4_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
lach4_tissues <- otu_over_time("Lachnospiraceae (OTU 4)" , otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_lach4_tissues.png", lach4_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)


#Examine "Lachnospiraceae (OTU 32)"  over time----
lach32_stools <- otu_over_time("Lachnospiraceae (OTU 32)" , otu_stools)+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/5_days_PEG_lach32_stools.png", lach32_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
lach32_tissues <- otu_over_time("Lachnospiraceae (OTU 32)" , otu_tissues)+
  scale_x_continuous(breaks = c(0, 4, 6, 20, 30),
                     limits = c(0,31),
                     minor_breaks = c(3.5, 4.5, 5.5, 6.5, 19.5, 20.5, 29.5, 30.5)) +
  theme(legend.position = "bottom")

save_plot(filename = "results/figures/5_days_PEG_lach32_tissues.png", lach32_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Customize days with scale_X_continuous



#Examine changes that happen in WMR group "5-day PEG 3350 + 10-day recovery" (baseline day -5 versus day 1)----
WMR_d0_d15_pairs <- otu_stools %>%
  filter(group == "WMR" & otu == "Bacteroides (OTU 1)") %>% #Limit to group "C" and randomly pick an OTU just to figure out what mice have sequence data
  filter(day == 1 | day == 15) %>%
  filter(duplicated(unique_mouse_id)) %>% #Pull mouse ids with sequence data for both day -1 and day 0
  pull(unique_mouse_id) #6 mice

#Wilcoxon signed rank test for all day 1, day 15 pairs at the OTU level:
otus_WMR_pairs <- otu_stools %>%
  filter(unique_mouse_id %in% WMR_d0_d15_pairs) %>% 
  filter(day == 1 | day == 15) %>% #Select timepoints to test
  group_by(otu) %>%
  nest() %>%
  mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>%
  mutate(median = map(data, get_rel_abund_median_day)) %>%
  unnest(c(model, median)) %>%
  ungroup()

#Adjust p-values for testing multiple OTUs
otus_WMR_pairs_stats_adjust <- otus_WMR_pairs %>%
  select(otu, statistic, p.value, method, alternative, `1`, `15`) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv(path = "data/process/5_days_PEG_otu_WMR_d1vd15.tsv")

#Create list of OTUs to plot for WMR group
otus_WMR_pairs_stats_adjust %>% filter(p.value.adj < 0.05) #No p-values survive multiple hypothesis correction
#Look at top 15 OTUs by unadjusted p-values
WMR_OTUs <- otus_WMR_pairs_stats_adjust %>% 
  arrange(p.value) %>% 
  slice_head(n=15) %>% 
  pull(otu)
#Add C. diff OTU to the list of OTUs to plot
WMR_OTUs <- c(WMR_OTUs, "Peptostreptococcaceae (OTU 12)")

#Heatmap of 5-day PEG + 10 day recovery (WMR) group----
#Create heatmaps of mock challenged mice----
#Mock stool samples
otu_WMR_stools  <-  otu_stools %>% 
  filter(group == "WMR")
hm_WMR_stool_days <- otu_WMR_stools %>% distinct(day) %>% pull(day)
hm_WMR_stool <-   otu_WMR_stools %>%
    mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
    mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                             "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the groups  
    filter(otu %in% WMR_OTUs) %>%
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
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 16))+ # Change font size for entire plot+
  scale_x_discrete(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30), labels = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_otus_heatmap_stools_WMR.png", hm_WMR_stool, base_height = 8, base_width = 8)


#Examine impacts of clindamycin and PEG3350 treatments on bacterial OTUs----
#Nov. 2020: Need to update this now that we have data from all timepoints

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

#Genus Level Analysis----
agg_genus_data_subset = five_day_PEG_subset(agg_genus_data) %>% 
  mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                           "5", "6", "7", "8", "9", "10", "15", "20", "25", "30"))

#Pairwise comparisons
#Create baseline WMR: Day -15 vs Day -10; Rest of groups Day -5 vs D1
# -15, -10, -5, -4 relevant timepoints for analysis of WMR 
WMR_baseline = agg_genus_data_subset %>%
  filter(group == "WMR", day == -15 | day == -10) %>%
  filter(duplicated(unique_mouse_id)) %>%
  pull(unique_mouse_id)
ro_groups_baseline = agg_genus_data_subset %>%
  filter(group != "WMR" | group != "C", day == -5 | day == 1) %>%
  filter(duplicated(unique_mouse_id)) %>%
  pull(unique_mouse_id)
#baseline_mice = rbind(WMR_baseline, ro_groups_baseline) %>%
 # filter(duplicated(unique_mouse_id)) %>%
 # pull(unique_mouse_id)

#Dataframes for statistical analysis at genus level
WMR_paired_genus <- agg_genus_data_subset %>%
  filter(unique_mouse_id %in% WMR_baseline) %>% 
  filter(day == -15 | day == -10) %>% 
  mutate(day = as.factor(day)) %>% 
  select(day, genus, agg_rel_abund)

ro_groups_paired_genus <- agg_genus_data_subset %>%
  filter(unique_mouse_id %in% ro_groups_baseline) %>% 
  filter(day == -5 | day == 1) %>% 
  mutate(day = as.factor(day)) %>% 
  select(day, genus, agg_rel_abund)

#Wilcoxon signed rank test for all pairs at the genus level: Must separate WMR from rest of groups due to different time points
WMR_genus_peg_effect_pairs <- WMR_paired_genus %>% 
  group_by(genus) %>% 
  nest() %>% 
  mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$day)) %>% tidy())) %>% 
  mutate(median = map(data, get_rel_abund_median_day)) %>% 
  unnest(c(model, median)) %>% 
  ungroup()

ro_groups_genus_peg_effect_pairs <- ro_groups_paired_genus %>% 
  group_by(genus) %>% 
  nest() %>% 
  mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$day)) %>% tidy())) %>% 
  mutate(median = map(data, get_rel_abund_median_day)) %>% 
  unnest(c(model, median)) %>% 
  ungroup()

#Adjust p-values for testing multiple genera
WMR_genus_effect_adj <- WMR_genus_peg_effect_pairs %>% 
  select(genus, statistic, p.value, method, `-15`, `-10`) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>%
  rename("baseline" = `-15`, "Post-PEG" = `-10`) #Rename to established baseline and post-PEG timepoints
ro_groups_effect_adj <- ro_groups_genus_peg_effect_pairs %>% 
  select(genus, statistic, p.value, method, `-5`, `1`) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>%
  rename("baseline" = `-5`, "Post-PEG" = `1`) #Rename to established baseline and post-PEG timepoints
  

#Pull significant genera from adjusted p-value dataframes
WMR_sig_genus_pairs <- pull_significant_taxa(WMR_genus_effect_adj, genus)
ro_groups_sig_genus_pairs <- pull_significant_taxa(ro_groups_effect_adj, genus)

#Combine and remove duplicated genera to find all significant genera before and after PEG
sig_genus_pairs <- c(WMR_sig_genus_pairs, ro_groups_sig_genus_pairs) %>%
  unique()

#Define Facet labels for faceting by genus in line plots
facet_labels <- sig_genus_pairs
names(facet_labels) <- sig_genus_pairs

#Find which days have samples from all groups
day_freq <- agg_genus_data_subset %>%
  group_by(day, group) %>%
  summarize(day = unique(day), group = color_groups) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq == 4)
  
genus_pair_days = unique(day_freq$day) #Days -1, -5, 0, 1, 3, 5, and 6 have samples from all groups

 #Plot all significant genus pairs for all timepoints and groups----
sig_genus_pair_plot <- agg_genus_data_subset %>% 
  filter(day %in% genus_pair_days) %>%
  mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
  mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                           "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% 
  filter(genus %in% sig_genus_pairs) %>%
  group_by(group, genus, day) %>%
  mutate(median_abund =median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%
  ggplot() +
  geom_line(mapping = aes(x = day, y = median_abund, group = group, color = group), alpha = 0.6, size = 1, show.legend = FALSE) +
  # geom_point(mapping = aes(x = day, y = agg_rel_abund, group = group, color = group), alpha = 0.6, size = 1, show.legend = FALSE) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels) +
  theme_classic()+
  labs(title=NULL,
       x="Days Post-Infection",
       y="Relative Abundance") +
  facet_wrap(~genus, labeller = labeller(genus = facet_labels)) +
  theme(plot.title=element_text(hjust=0.5),
         strip.background = element_blank(), #get rid of box around facet_wrap labels
         axis.text.y = element_markdown(), #Have only the OTU names show up as italics
         text = element_text(size = 16)) # Change font size for entire plot
save_plot(filename = "results/figures/5_days_PEG_genus_pairs.png", sig_genus_pair_plot, base_height = 7, base_width = 20)





#Plot of significant genera post-PEG treatment----
facet_labels <- c("Baseline", "PEG Treatment")
names(facet_labels) <- c(-5, 1)

#Plot all groups but WMR group (has different baseline and post-treatment timepoints)
peg_impacted_genera_no_WMR_plot <- agg_genus_data_subset %>% 
  filter(genus %in% ro_groups_sig_genus_pairs) %>% 
  filter(day %in% c(-5, 1)) %>% #Select baseline and post-peg timepoints
  filter(group != "WMR") %>% # Remove WMR as these timepoints do not represent its baseline and post-treatment timepoints
  mutate(genus=factor(genus, levels=ro_groups_sig_genus_pairs)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/5437) %>% 
  ggplot(aes(x= genus, y=agg_rel_abund, color=group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  geom_hline(yintercept=1/5437, color="gray")+
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
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  facet_wrap(~day, labeller = labeller(day = facet_labels), scales = "fixed")+
  theme(plot.title=element_text(hjust=0.5),
        text = element_text(size = 16),# Change font size for entire plot
        axis.text.y = element_markdown(face = "italic"), #Make sure genera names are in italics
        strip.background = element_blank(),
        legend.position = "none") 
save_plot(filename = paste0("results/figures/5_days_peg_impacted_genera_no_WMR_plot.png"), peg_impacted_genera_no_WMR_plot, base_height = 9, base_width = 9)

facet_labels <- c("Baseline", "PEG Treatment")
names(facet_labels) <- c(-15, -10)

#Plot PEG's effect on the WMR group alone
peg_impacted_genera_WMR_plot <- agg_genus_data_subset %>%
  filter(genus %in% WMR_sig_genus_pairs) %>% 
  filter(day %in% c(-15, -10)) %>% #Select baseline and post-peg timepoints
  mutate(genus=factor(genus, levels=WMR_sig_genus_pairs)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/5437) %>% 
  ggplot(aes(x= genus, y=agg_rel_abund, color=group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  geom_hline(yintercept=1/5437, color="gray")+
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
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  facet_wrap(~day, labeller = labeller(day = facet_labels), scales = "fixed")+
  theme(plot.title=element_text(hjust=0.5),
        text = element_text(size = 16),# Change font size for entire plot
        axis.text.y = element_markdown(face = "italic"), #Make sure genera names are in italics
        strip.background = element_blank(),
        legend.position = "none") 
save_plot(filename = paste0("results/figures/5_days_peg_impacted_genera_WMR_plot.png"), peg_impacted_genera_WMR_plot, base_height = 9, base_width = 9)

  



