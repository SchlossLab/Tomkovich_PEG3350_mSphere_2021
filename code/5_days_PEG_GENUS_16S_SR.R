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

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Remove environmental variables that are no longer needed
rm(agg_genus, agg_taxa_data, duplicated_seq_samples, duplicates_to_drop, genus_data, peg3350.files,
   contaminated_notes, contaminated_samples, seq_files_missing_from_metadata, prep_notes)

#Genus level analysis----
#Subset genus data to just the 5-day PEG subset and separate by sample type (stools versus tissues)
genus_subset <- five_day_PEG_subset(agg_genus_data)
view(genus_subset)
#Create subset dataframes of the 5-days PEG diversity data for just stool samples, tissues. 
genus_stools <- subset_stool(genus_subset)
genus_tissues <- subset_tissue(genus_subset)
#The above subsets exclude mock challenged mice (group = WMN or CN)
#Also create dataframes of diversity data that includes mock challenged mice (WMN and C), separated into stool and tissue samples
genus_mock_stools <- subset_stool(add_mocks(genus_subset, agg_genus_data))
genus_mock_tissues <- subset_tissue(add_mocks(genus_subset, agg_genus_data))

#Experimental days to analyze with the Kruskal-Wallis test (timepoints with 16S data for at least 3 groups)
#Baseline (before treatment) for WMR is day -15. For C, WM, and WMC baseline is day -5
stool_test_days <- c(-5, -1, 0, 1, 2, 3, 4, 5, 6, 10, 30)

tissue_test_days <- c(6, 30) #Only 2 days with samples from at least 3 groups

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



#Genus at varied across treatment groups and were shared across days 
#Stools
shared_sig_stools_genus_D1toD6 <- intersect_all(sig_genus_stools_day1, sig_genus_stools_day2,                
                                                sig_genus_stools_day3, sig_genus_stools_day4, 
                                               sig_genus_stools_day5, sig_genus_tissues_day6) #fill in different days to compare
view(shared_sig_stools_genus_D1toD6)
print(shared_sig_stools_genus_D1toD6)

# "Bacteroides"                        "Peptostreptococcaceae Unclassified"
# "Enterobacteriaceae Unclassified"    "Ruminococcaceae Unclassified"


#Tissues
shared_sig_tissues_genus <- intersect_all(sig_genus_tissues_day30, sig_genus_tissues_day6) #fill in different days to compare
print(shared_sig_tissues_genus)

# "Unclassified"                 "Clostridiales Unclassified"   "Firmicutes Unclassified"     
# "Clostridium XlVb"             "Bacteroides"                  "Ruminococcaceae Unclassified"
# "Butyricicoccus"

#Function to plot a list of genera across groups of mice at a specific timepoint:
#Arguments: 
#sample_df = subset dataframe of samples to be plotted
#genus = list of genera to plot
#timepoint = day of the experiment to plot
plot_genus_dx <- function(sample_df, genus, timepoint){
  sample_df %>%
    filter(genus %in% genus) %>%
    filter(day == timepoint) %>%
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% # 2,000 is 2 times the subsampling parameter of 1000
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
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "none",
          axis.text.y = element_markdown(), #Have only the genus names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Plots of the relative abundances of Genera that significantly varied across sources of mice from day -1 to day 1----
genus_d1 <- plot_genus_dx(genus_stools, sig_genus_stools_day1 [1:20], 1)+
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Day 1 post-infection Stools")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d1_top20.png", genus_d1, base_height = 7, base_width = 8)

genus_d4 <- plot_genus_dx(genus_stools, sig_genus_stools_day4 [1:20], 4)+
  ggtitle("Day 4 post-infection Stools")+ #Title plot
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d4_top20.png", genus_d4, base_height = 7, base_width = 8)

genus_d10 <- plot_genus_dx(genus_stools, sig_genus_stools_day10[1:20], 10)+
  ggtitle("Day 10 post-infection Stools")+ #Title plot
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d10.png", genus_d10, base_height = 7, base_width = 8)

genus_d6 <- plot_genus_dx(genus_stools, `sig_genus_stools_day6` [1:20], 6)+
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  ggtitle("Day 6 post-infection Stools")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_stools_d6_top20.png", genus_d6, base_height = 7, base_width = 8)

genus_tissues_d6 <- plot_genus_dx(genus_tissues, sig_genus_tissues_day6[1:20], 6)+
  ggtitle("Day 6 post-infection Tissue")+ #Title plot
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_tissues_d6.png", genus_tissues_d6, base_height = 7, base_width = 8)

genus_tissues_d30 <- plot_genus_dx(genus_tissues, sig_genus_tissues_day30[1:20], 30)+
  ggtitle("Day 30 post-infection Tissue")+ #Title plot
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate Genera
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/5_days_PEG_genus_tissues_d30.png", genus_tissues_d30, base_height = 7, base_width = 8)

#Heatmaps of significant Genera over time facet by group----
#Create list of Genera to plot
all_sig_genus <- c(`sig_genus_stools_day-5`, `sig_genus_stools_day-1`, sig_genus_stools_day0, sig_genus_stools_day1, 
                  `sig_genus_stools_day2`, `sig_genus_stools_day3`, sig_genus_stools_day4, sig_genus_stools_day5,
                  `sig_genus_stools_day6`, `sig_genus_stools_day10`, sig_genus_stools_day30)
#370 total Genera
unique_sig_genus <-unique(all_sig_genus)
#130 unique significant Genera
#Rank the 130 signifiant Genera in order of agg_rel_abund, select the top 20
hm_sig_genus_abund <-  genus_stools %>% 
  filter(genus %in% unique_sig_genus) %>%
  group_by(genus) %>% 
  summarize(median=(median(agg_rel_abund + 1/2000))) %>% #Median relative abundance for all samples for a particular genus
  arrange(desc(median)) %>% #Arrange largest to smallest
  slice_head(n = 20) %>% 
  pull(genus)

#Rank Genera by adjusted p-value
hm_sig_genus_p_adj <- kw_genus_stools %>% 
  filter(p.value.adj < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(genus) %>% 
  slice_head(n = 25) %>% 
  pull(genus)

hm_stool_days <- diversity_stools %>% distinct(day) %>% pull(day)
facet_labels <- color_labels #Create descriptive labels for facets
names(facet_labels) <- c("C", "WM", "WMC", "WMR") #values that correspond to group, which is the variable we're faceting by
hm_stool <- hm_plot_genus(genus_stools, hm_sig_genus_p_adj, hm_stool_days)+
  scale_x_discrete(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30), labels = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_genus_heatmap_stools.png", hm_stool, base_height = 14, base_width = 15)

#Create heatmap of significant Genera for tissue samples----
#Rank Genera by adjusted p-value
hm_sig_genus_p_adj_tissues <- kw_genus_tissues %>% 
  filter(p.value.adj < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(genus) %>% 
  slice_head(n = 25) %>% 
  pull(genus)
hm_tissue_days <- diversity_tissues %>% distinct(day) %>% pull(day)
hm_tissues <- hm_plot_genus(genus_tissues, hm_sig_genus_p_adj_tissues, hm_tissue_days)+
  scale_x_discrete(breaks = c(4, 6, 20, 30), labels = c(4, 6, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_genus_heatmap_tissues.png", hm_tissues, base_height = 14, base_width = 15)
#Examine Genera that were significant in stool samples
hm_tissues_stool_genus <- hm_plot_genus(genus_tissues, hm_sig_genus_p_adj, hm_tissue_days)+
  scale_x_discrete(breaks = c(4, 6, 20, 30), labels = c(4, 6, 20, 30)) 
save_plot(filename = "results/figures/5_days_PEG_genus_heatmap_tissues_stool_genus.png", hm_tissues_stool_genus, base_height = 14, base_width = 15)
#Pull list of Genera that overlap between stool & tissue heatmaps
hm_overlap <- intersect_all(hm_sig_genus_p_adj, hm_sig_genus_p_adj_tissues)
#11 Genera overlap:  
#"Lachnospiraceae (OTU 33)"       "Blautia (OTU 19)"              
#"Ruminococcaceae (OTU 50)"       "Ruminococcaceae (OTU 54)"      
#"Ruminococcaceae (OTU 92)"       "Oscillibacter (OTU 45)"        
# "Lachnospiraceae (OTU 30)"       "Lachnospiraceae (OTU 31)"      
# "Lachnospiraceae (OTU 4)"        "Peptostreptococcaceae (OTU 12)"
# "Enterobacteriaceae (OTU 2)" 

