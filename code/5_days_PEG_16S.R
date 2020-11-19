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

metadata <- metadata %>%
  mutate(day = as.integer(day))  #Day variable (transformed to integer to get rid of decimals on PCoA animation

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Create functions to further subset our data by sample type and or add mock challenged mice (CN and WMN groups)----
#Function to filter dataframe (df) to examine stool samples
subset_stool <- function(df){
  df %>% filter(sample_type == "stool")
}
#Function to filter dataframe (df) to examine tissue samples
subset_tissue <- function(df){
  df %>% filter(!sample_type == "stool")
}

#Function to add the mock challenged mice (CN and WMN groups) to a dataframe (diversity_data or agg_otu_data)
add_mocks <- function(subset_df, original_df){
  subset_df %>% 
  add_row(original_df %>% filter(group %in% c("CN", "WMN"))) #Add the groups of mock challenged mice
}

#Function to figure out how many samples we have per group per day for each subset (subset = dataframe of a subset of samples from the 5-days PEG subset)
count_subset <- function(subset){
  subset %>% 
    group_by(group) %>%
    count(day) %>% 
    arrange(day)
}

#Alpha diversity analysis----

#Subset diversity data to just the 5-day PEG subset:
diversity_subset <- five_day_PEG_subset(diversity_data)
#five_day_PEG_subset() will exclude mock challenged mice (group = WMN or CN)

#Create subset dataframes of the 5-days PEG diversity data for just stool samples, tissues. 
diversity_stools <- subset_stool(diversity_subset)
diversity_tissues <- subset_tissue(diversity_subset)
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
stool_mock_test_days
tissue_test_days <- c(6, 30) #Only 2 days with samples from at least 3 groups
tissue_mock_test_days

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
sig_richness_days_stissues <- pull_sig_days(kw_richness_tissues)

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

save_plot(filename = paste0("results/figures/5_Day_PEG_PCoA.png"), pcoa_subset_plot, base_height = 5, base_width = 4.5)

#Pull 5_Day_PEG subset of PCoA 
#Tissue subset
pcoa_5_day_PEG_tissues <- read_tsv("data/process/5_day_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and use left_join to keep all samples in pcoa data frame
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

pcoa_5_day_PEG_tissues<- subset(pcoa_5_day_PEG_tissues, group %in% c("C", "WM", "WMN", "WMC"))

#Pull axes from loadings file
pcoa_axes_5_day_PEG_tissues <- read_tsv("data/process/5_day_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG_tissues %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG_tissues %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot_tissue <- plot_pcoa(pcoa_5_day_PEG_tissues)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))

save_plot(filename = paste0("results/figures/5_Day_PEG_tissues_PCoA.png"), pcoa_subset_plot_tissue, base_height = 5, base_width = 4.5)

#Stool Subset
pcoa_5_day_PEG_stool <- read_tsv("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

pcoa_5_day_PEG_stool<- subset(pcoa_5_day_PEG_stool, group %in% c("C", "WM", "WMN", "WMC"))

pcoa_axes_5_day_PEG_stool <- read_tsv("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_5_day_PEG_stool %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot_stool <- plot_pcoa(pcoa_5_day_PEG_stool)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))

save_plot(filename = paste0("results/figures/5_Day_PEG_stool_PCoA.png"), pcoa_subset_plot_stool, base_height = 5, base_width = 4.5)

#Remove legend
pcoa_plot_time <- plot_pcoa(pcoa_5_day_PEG_tissues)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))+
  theme(legend.position = "none")+ #remove legend
  facet_wrap(~ day)

#PCoAs of select timepoints of interst

#Animation of PCoA plot over time for all sequenced samples ----
#Source: Will Close's Code Club from 4/12/2020 on plot animation
pcoa_animated <- plot_pcoa(pcoa_5_day_PEG_tissues)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))+
  labs(title = 'Day: {frame_time}') + #Adds time variable to title
  transition_time(day)+  #Day variable used to cycle through time on animation
  shadow_mark() #Shows previous timepoints

# Implement better frames per second for animation
pcoa_gif <- animate(pcoa_animated, duration = 6, fps = 10,
                    res = 150, width = 20, height = 20, unit = "cm")

# Save as gif file
anim_save(animation = pcoa_gif, filename = 'results/5_days_PEG_pcoa_over_time_tissues.gif')


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

# Perform kruskal wallis tests at the otu level for the stool samples----
for (d in stool_test_days){
  kruskal_wallis_otu(d, otu_stools, "stools")
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/5_days_PEG_otu_stats_day_", d, "_stools.tsv"))
  name <- paste("sig_otu_stools_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, otu))
}

# Perform kruskal wallis tests at the otu level for the tissue samples----
for (d in tissue_test_days){
  kruskal_wallis_otu(d, otu_tissues, "tissues")
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/5_days_PEG_otu_stats_day_", d, "_tissues.tsv"))
  name <- paste("sig_otu_tissues_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, otu))
}

#OTUs that varied across treatment groups and were shared across days 
shared_sig_otus <- intersect_all() #fill in different days to compare

#Function to plot a list of OTUs across sources of mice at a specific timepoint:
#Arguments: otus = list of otus to plot; timepoint = day of the experiment to plot
plot_otus_dx <- function(otus, timepoint){
  agg_otu_data %>%
    filter(otu %in% otus) %>%
    filter(day == timepoint) %>%
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% # 2,000 is 2 times the subsampling parameter of 1000
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

#Function to plot an otu_over_time
#otu_plot = otu to plot in quotes. Ex: "Peptostreptococcaceae (OTU 12)"
#sample_df = subset dataframe of just stool or tissue samples
otu_over_time <- function(otu_plot, sample_df){
  specify_otu_name <- sample_df %>% 
    filter(otu == otu_plot) %>% 
    pull(otu_name)
  otu_median <- sample_df %>% 
    filter(otu == otu_plot) %>% 
    group_by(group, day) %>% 
    summarize(median=(median(agg_rel_abund + 1/2000))) %>% 
    ungroup
  otu_mice <- sample_df %>% 
    filter(otu == otu_plot) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>%
    select(day, agg_rel_abund, otu, group)
  otu_time <- ggplot(NULL)+
    geom_point(otu_mice, mapping = aes(x=day, y=agg_rel_abund, color=group), size  = 1.5, position = position_dodge(width = 0.6))+
    geom_line(otu_median, mapping = aes(x=day, y=median, color=group), size = 1, show.legend = FALSE)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    geom_hline(yintercept=1/1000, color="gray")+
    labs(title=specify_otu_name,
         x="Day",
         y="Relative abundance (%)") +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_markdown(hjust = 0.5),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"),  # Add gray lines to clearly separate symbols by days)
          text = element_text(size = 18)) # Change font size for entire plot
}
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


#Customize days with scale_X_continuous

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
