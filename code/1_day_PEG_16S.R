source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match 1 Day Peg Plots
color_scheme <- c("#225ea8", "#238b45", "#88419d") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c( "1RM1", "C", "M1")
color_labels <- c( "1-day PEG 3350 + 1-day recovery", "Clind.", "1-day PEG 3350")

#Subset alpha diversity data (16S_common_files) to analyze one day PEG subset mice----
diversity_data_subset <- one_day_PEG_subset(diversity_data) %>%
  mutate(day = replace(day, day == -1, "PT"),
         day = replace(day, day == -2, "PT")) #Replace day -2 and day -1 with PT to represent the pretreatment timepoint

 diversity_data_subset <- diversity_data_subset %>%
   filter(day %in% c("PT", 0, 1, 2, 4, 5, 7)) %>%
   mutate(day = fct_relevel(day, "PT", "0", "1" , "2" , "4", "5" , "7")) #Specify the order of the days for plotting overtime
 
#Get exp days sequenced in both pretreatment and other days
exp_days_seq <- unique(diversity_data_subset %>% pull(day)) 

#Remove Days 0 and 4 because not all three groups are represented
exp_days_seq <- exp_days_seq[1:5]  #Removes last two positions of day 0 and day 4

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Function to test at the otu level:
agg_otu_data_subset <- one_day_PEG_subset(agg_otu_data) %>%
  mutate(day = replace(day, day == -1, "PT"),
         day = replace(day, day == -2, "PT")) %>% #Replace day -2 and day -1 with PT to represent the pretreatment timepoint
  mutate(day = fct_relevel(day, "PT", "0", "1" , "2" , "4", "5" , "7")) #Specify the order of the days for plotting overtime

#Kruskal-Wallis Function to test for differences in Shannon Diveristy Index across groups with Benjamini Hochberg correction----
#Adapted from 5_days_PEG_16S.R script, removed subset_name argument as only one sample type in 1 Day subset
#Arguments: 
#diversity_subset <- subset (stools or tissue samples) of diversity_data to perform statistical test on
#timepoint = timepoints to assess differences between groups specific to the subset (stool or tissue)
kruskal_wallis_shannon <- function(diversity_subset, timepoints){
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
    select(day, statistic, p.value, parameter, method, "C", "1RM1", "M1") %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/1_days_PEG_shannon_stats_subset.tsv"))
}
#Test with shannon (Only stool samples in 1 Day Subset)
kw_shannon <- kruskal_wallis_shannon(diversity_data_subset, exp_days_seq)
sig_shannon_days <- pull_sig_days(kw_shannon)

#Kruskal-Wallis Function to test for differences in richness (sobs) between groups on a particular day with Benjamini Hochberg correction----
#Adapted from 5_days_PEG_16S.R script, removed subset_name argument as only one sample type in 1 Day subset
#Arguments: 
#diversity_subset <- subset (stools or tissue samples) of diversity_data to perform statistical test on
#timepoint = timepoints to assess differences between groups specific to the subset (stool or tissue)
kruskal_wallis_richness <- function(diversity_subset, timepoints){
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
    select(day, statistic, p.value, parameter, method, "C", "1RM1", "M1") %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/1_days_PEG_richness_stats_subset.tsv"))
}
#Test with richness (Only stool samples in 1 Day subset)
kw_richness <- kruskal_wallis_richness(diversity_data_subset, exp_days_seq)
sig_richness_days <- pull_sig_days(kw_richness) 

#Shannon and richness plots for the 1 Day PEG subset----
#Statistical annotation labels:
x_annotation <- sig_shannon_days
y_position <- max(diversity_data_subset$shannon)+0.15
label <- kw_label(kw_shannon)
#Plot
shannon_stools <- plot_shannon_overtime(diversity_data_subset) +
  scale_x_discrete(breaks = c("PT", 0, 1, 2, 4, 5, 7)) +
  geom_vline(xintercept = c((1:7) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/1_Day_PEG_shannon.png", shannon_stools, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Statistical annotation labels:
x_annotation <- sig_richness_days
y_position <- max(diversity_data_subset$sobs)+5
label <- kw_label(kw_richness)
#Plot
richness_tissues <- plot_richness_overtime(diversity_data_subset) +
  scale_x_discrete(breaks = c("PT", 0, 1, 2, 4, 5, 7)) +
  geom_vline(xintercept = c((1:7) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/1_Day_PEG_richness.png", richness_tissues, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#PCoa Analysis----
#Pull 1_Day_PEG subset of PCoA data
pcoa_1_day_PEG <- read_tsv("data/process/1_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames 
  mutate(day = as.integer(day)) %>% #Day variable (transformed to integer to get rid of decimals on PCoA animation
  filter(!is.na(axis1)) %>% #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff
  filter(day > -3)#limit to experimental time frameampling cutoff

#Pull axes from loadings file
pcoa_axes_1_day_PEG <- read_tsv("data/process/1_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_1_day_PEG %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_1_day_PEG %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa_subset_plot <- plot_pcoa(pcoa_1_day_PEG)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))
save_plot(filename = paste0("results/figures/1_Day_PEG_PCoA.png"), pcoa_subset_plot, base_height = 5, base_width = 4.5)

#Create stand alone legend
group_legend <- pcoa_1_day_PEG  %>%
  ggplot(aes(x = axis1, y = axis2, color = group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  geom_point()+ theme_classic()
group_legend <- get_legend(group_legend)
save_plot("results/figures/1_day_PEG_pcoa_legend.png", group_legend, base_height = .8, base_width = 2.2)

pcoa_animated <- plot_pcoa(pcoa_1_day_PEG)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  labs(title = 'Day: {frame_time}') + #Adds time variable to title
  transition_time(day)+  #Day variable used to cycle through time on animation
  shadow_mark() #Shows previous timepoints

# Implement better frames per second for animation
pcoa_gif <- animate(pcoa_animated, duration = 7, fps = 10,
                    res = 150, width = 20, height = 20, unit = "cm")

# Save as gif file
anim_save(animation = pcoa_gif, filename = 'results/1_Day_PEG_pcoa_over_time.gif')


#Kruskal-Wallis Function to test for differences in OTUs across groups on a particular day----
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
    select(otu, statistic, p.value, parameter, method, "C", "1RM1", "M1") %>%
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
    arrange(p.value.adj) %>%
    write_tsv(path = paste0("data/process/1_Day_PEG_otu_stats_day_", timepoint, ".tsv"))
}

#Create empty data frame to combine stat dataframes for all days that were tested
kw_otu <- data.frame(otu=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                            C =double(), `1RM1` =double(), M1=double(),
                            p.value.adj=double(),day=character()) %>% 
  rename(`1RM1` = X1RM1) #Rename 1RM1 to get rid of X
# Perform kruskal wallis tests at the otu level for all days of the experiment that were sequenced
for (d in exp_days_seq){
  kruskal_wallis_otu(d)
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/1_Day_PEG_otu_stats_day_", d, ".tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_otu_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, otu))
  kw_otu <- add_row(kw_otu, stats)  #combine all the dataframes together
}

#Shared significant genera across from Days 1, 2 and 7----
shared_sig_otus_D1toD7 <- intersect_all(`sig_otu_day1`, sig_otu_day2, sig_otu_day5, sig_otu_day7, sig_otu_dayPT)
view(shared_sig_otus_D1toD7) #0
print(shared_sig_otus_D1toD7)
#0 OTUs overlap

#No taxa significant 
for (d in exp_days_seq){
  #Make a list of top 10 otus across sources of mice for a specific day
  top10 <- read_tsv(file = paste0("data/process/1_Day_PEG_otu_stats_day_", d, ".tsv")) %>% 
    slice_head(n = 10) %>%  #Select the top 10 otus (already arranged by adjusted p values)
    pull(otu)
  name <- paste("top10_otu_day", d, sep = "")
  assign(name, top10)
}

#Plot top 10 taxa for each day (no significant differences in the relative abundances between groups:
#Function to plot a list of OTUs across sources of mice at a specific timepoint:
#Arguments: otus = list of otus to plot; timepoint = day of the experiment to plot
plot_otus_dx <- function(otus, timepoint){
  agg_otu_data_subset %>% 
    filter(otu %in% otus) %>%
    filter(day == timepoint) %>%
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% # 2000 is 2 times the subsampling parameter of 1000
    ggplot(aes(x= otu_name, y=agg_rel_abund, color=group))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    geom_hline(yintercept=1/1000, color="gray")+ #Limit of detection = 1/our subsampling parameter
    stat_summary(fun = 'median',
                 fun.max = function(x) quantile(x, 0.75),
                 fun.min = function(x) quantile(x, 0.25),
                 position = position_dodge(width = 1)) +
    labs(title=NULL,
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/2000, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "bottom",
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Plots of the top 20 OTUs that varied across sources at each timepoint----
#Pre-treatment (PT): top 10 OTUs, no difference between groups
PT_otus <- plot_otus_dx(top10_otu_dayPT, "PT")+ #Pick top 10 OTUs
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_PT_ns_otus.png", PT_otus, base_height = 9, base_width = 7)
#Day 1: top 10 OTUs, no difference between groups
D1_otus <- plot_otus_dx(top10_otu_day1, 1) +#Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_D1_ns_otus.png", D1_otus, base_height = 9, base_width = 7)
#Day 2: top 10 OTUs, no difference between groups
D2_otus <- plot_otus_dx(top10_otu_day2, 2) + #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_D2_ns_otus.png", D2_otus, base_height = 9, base_width = 7)
#Day 5: top 10 OTUs, no difference between groups
D5_otus <- plot_otus_dx(top10_otu_day5, 5)+ #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_D5_ns_otus.png", D5_otus, base_height = 9, base_width = 7)
#Day 7: top 10 OTUs, no difference between groups
D7_otus <- plot_otus_dx(top10_otu_day7, 7)+ #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_D7_ns_otus.png", D7_otus, base_height = 9, base_width = 7)

#Examine C. difficile otu over time----
peptostrep <- otu_over_time("Peptostreptococcaceae (OTU 12)", agg_otu_data_subset)+
  scale_x_discrete(breaks = c("PT", 0, 1, 2, 4, 5, 7)) +
  geom_vline(xintercept = c((1:7) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "bottom")
save_plot(filename = "results/figures/1_day_PEG_otu_peptostreptococcaceae.png", peptostrep, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Heatmaps of OTUs (lowest P value since none were significant) over time, facet by group----
#Rank OTUs by adjusted p-value
hm_otus_p <- kw_otu %>% 
  filter(p.value < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(otu) %>% 
  slice_head(n = 20) %>% 
  pull(otu)

hm_days <- c("PT", "0", "1", "2", "4", "5", "7")
facet_labels <- color_labels #Create descriptive labels for facets
names(facet_labels) <- c("1RM1", "C", "M1") #values that correspond to group, which is the variable we're faceting by
hm_otu <- hm_plot_otus(agg_otu_data_subset, hm_otus_p, hm_days)+
  scale_x_discrete(limits = c("PT", "0", "1", "2", "4", "5", "7"), breaks = c("PT", "0", "1", "2", "4", "5", "7"), labels = c("PT", "0", "1", "2", "4", "5", "7")) 
save_plot(filename = "results/figures/1_day_PEG_otus_heatmap.png", hm_otu, base_height = 10, base_width = 15)

#Plot heat map of OTUs that were significant in the 5 day subset
hm_5_day_otus <- c("Phenylobacterium (OTU 332)", "Lachnospiraceae (OTU 33)", "Ruminococcaceae (OTU 37)",
                   "Ruminococcaceae (OTU 98)", "Ruminococcaceae (OTU 65)", "Clostridium (OTU 51)",
                   "Bacteroides (OTU 1)", "Lactobacillus (OTU 13)", "Blautia (OTU 19)", 
                   "Ruminococcaceae (OTU 50)", "Ruminococcaceae (OTU 54)", "Ruminococcaceae (OTU 92)",
                   "Bifidobacterium (OTU 28)", "Oscillibacter (OTU 45)", "Lachnospiraceae (OTU 30)", 
                   "Lachnospiraceae (OTU 29)", "Lachnospiraceae (OTU 31)", "Lactobacillus (OTU 23)", 
                   "Lachnospiraceae (OTU 4)", "Peptostreptococcaceae (OTU 12)", "Lachnospiraceae (OTU 32)", 
                   "Enterobacteriaceae (OTU 2)", "Clostridium (OTU 22)", "Lactobacillus (OTU 9)", "Lachnospiraceae (OTU 16)")
hm_5_day_otus_plot <- hm_plot_otus(agg_otu_data_subset, hm_5_day_otus, hm_days)+
  scale_x_discrete(limits = c("PT", "0", "1", "2", "4", "5", "7"), breaks = c("PT", "0", "1", "2", "4", "5", "7"), labels = c("PT", "0", "1", "2", "4", "5", "7")) 
save_plot(filename = "results/figures/1_day_PEG_otus_heatmap_5_day_otus.png", hm_5_day_otus_plot, base_height = 10, base_width = 15)

