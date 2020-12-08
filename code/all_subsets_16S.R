source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

metadata <- metadata %>%
  mutate(day = as.integer(day))  #Day variable (transformed to integer to get rid of decimals on PCoA animation
set.seed(19760620) #Same seed used for mothur analysis

#Compare the stool and tissue samples of all mice that were infected with C. difficile----
#Define groups:
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8", ) #need to modify if we decide to show all groups together
all_groups <- c("C", "WM", "WMC", "WMR",
                "M1", "1RM1",
                "CWM", "FRM", "RM", 
                "CN", "WMN")
all_labels <- c( "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 + 10-day recovery",
                 "1-day PEG 3350", "1-day PEG 3350 + 1-day recovery",
                 "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350",
                 "Clind. without infection", "5-day PEG 3350 without infection")

#Figure out days to test in stool & tissue subsets----
all_diversity <- diversity_data %>% 
  filter(group %in% all_groups)
all_diversity_stools <- subset_stool(all_diversity)
all_diversity_tissues <- subset_tissue(all_diversity)
num_stool <- count_subset(all_diversity_stools)
num_tissue <- count_subset(all_diversity_tissues)
stool_test_days <- c(-15, -5, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 30)
tissue_test_days <- c(4, 6, 30) #Only days with samples from at least 3 groups

#Test for differences in OTU relative abundances across groups at specific timepoiints----
#subset agg_otu_data
all_otu <- agg_otu_data %>% 
  filter(group %in% all_groups)
all_otu_stools <- subset_stool(all_otu)
all_otu_tissues <- subset_tissue(all_otu)
#Function to test stool or tissue samples at the otu level:
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
    write_tsv(path = paste0("data/process/all_otu_stats_day_", timepoint, "_", sample_type, ".tsv"))
}

#Create empty data frame to combine stat dataframes for all days that were tested
kw_otu_stools <- data.frame(otu=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                            C =double(), WM =double(),WMC=double(),WMR =double(),M1 =double(),`1RM1` =double(),
                            CWM =double(),FRM=double(),RM =double(),
                            CN=double(),WMN =double(),
                            p.value.adj=double(),day=double())%>% 
  rename(`1RM1` = X1RM1) #Rename 1RM1 to get rid of X

# Perform kruskal wallis tests at the otu level for all mice stool samples----
for (d in stool_test_days){
  kruskal_wallis_otu(d, all_otu_stools, "stools")
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/all_otu_stats_day_", d, "_stools.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_otu_stools_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, otu))
  kw_otu_stools <- add_row(kw_otu_stools, stats)  #combine all the dataframes together
}

#Create empty data frame to combine stat dataframes for all days that were tested
kw_otu_tissues <- data.frame(otu=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                             C =double(), WM =double(),WMC=double(),WMR =double(),
                             CWM =double(),FRM=double(),RM =double(),
                             CN=double(),WMN =double(),
                             p.value.adj=double(),day=double())

# Perform kruskal wallis tests at the otu level for the tissue samples----
for (d in tissue_test_days){
  kruskal_wallis_otu(d, all_otu_tissues, "tissues")
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/all_otu_stats_day_", d, "_tissues.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_otu_tissues_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, otu))
  kw_otu_tissues <- add_row(kw_otu_tissues, stats)  #combine all the dataframes together
}

#Create heatmap of significant OTUs for all stool samples----
#Rank OTUs by adjusted p-value
hm_sig_otus_p_adj <- kw_otu_stools %>% 
  filter(p.value.adj < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(otu) %>% 
  slice_head(n = 25) %>% 
  pull(otu)

hm_stool_days <- all_diversity_stools %>% distinct(day) %>% pull(day)
facet_labels <- all_labels #Create descriptive labels for facets
names(facet_labels) <- c("C", "WM", "WMC", "WMR",
                         "M1", "1RM1",
                         "CWM", "FRM", "RM",
                         "CN", "WMN") #values that correspond to group, which is the variable we're faceting by
hm_stool <- hm_plot_otus(all_otu_stools, hm_sig_otus_p_adj, hm_stool_days)+
  scale_x_discrete(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30), labels = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30)) 
save_plot(filename = "results/figures/all_otus_heatmap_stools.png", hm_stool, base_height = 15, base_width = 18)

#Create heatmap of significant OTUs for all tissue samples----
#Rank OTUs by adjusted p-value
hm_sig_otus_p_adj_tissues <- kw_otu_tissues %>% 
  filter(p.value.adj < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(otu) %>% 
  slice_head(n = 25) %>% 
  pull(otu)
hm_tissue_days <- all_diversity_tissues %>% distinct(day) %>% pull(day)
hm_tissues <- hm_plot_otus(all_otu_tissues, hm_sig_otus_p_adj_tissues, hm_tissue_days)+
  scale_x_discrete(breaks = c(4, 6, 20, 30), labels = c(4, 6, 20, 30)) 
save_plot(filename = "results/figures/all_otus_heatmap_tissues.png", hm_tissues, base_height = 14, base_width = 15)
#Pull list of otus that overlap between stool & tissue heatmaps
hm_overlap <- intersect_all(hm_sig_otus_p_adj, hm_sig_otus_p_adj_tissues)
# "Enterobacteriaceae (OTU 2)"     "Lachnospiraceae (OTU 33)"      
# "Peptostreptococcaceae (OTU 12)" "Blautia (OTU 19)"              
# "Ruminococcaceae (OTU 27)"       "Lachnospiraceae (OTU 31)"  