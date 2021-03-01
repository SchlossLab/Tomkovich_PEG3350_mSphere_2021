source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files
library(plotly)
library(viridis)

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

#Test for differences in OTU relative abundances across groups at specific timepoints----
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

# Perform pairwise Wilcoxan rank sum tests for otus that were significantly different across sources of mice on a series of days----
pairwise_day_otu <- function(sample_df, sample_type, timepoint, sig_otu_dayX){
  otu_stats <- sample_df %>% 
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
    write_tsv(path = paste0("data/process/all_otu_stats_day_", timepoint, "_", sample_type, "_sig.tsv"))
  #Format pairwise stats to use with ggpubr package
  plot_format_stats <- pairwise_stats %>% 
    select(-method) %>% 
    group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
    lapply(tidy_pairwise_otu) %>% 
    bind_rows() %>% 
    filter(!is.na(group2)) %>% #Remove rows that don't have a 2nd comparison group
    arrange(p.adj) %>% #Arrange by adjusted p value column 
    mutate(day = timepoint)
  return(plot_format_stats)  
}

#Pairwise Tests of OTUs that were significant on each day tested
#List of significant days:
pairwise_days_stools <- kw_otu_stools %>% filter(p.value.adj < 0.05) %>% 
  distinct(day) %>% pull(day)
#Pairwise test of significant days and their corresponding OTUs for stool samples 
#Create empty placeholder data frame
pairwise_otu_stools <- data.frame(otu=character(), group1=character(), group2=character(),
                                  p.adj = double(), day= double())
otu_dayn5_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", -5, `sig_otu_stools_day-5`)
otu_dayn1_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", -1, `sig_otu_stools_day-1`)
otu_day0_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 0, sig_otu_stools_day0)
otu_day1_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 1, sig_otu_stools_day1)
otu_day2_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 2, sig_otu_stools_day2)
otu_day3_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 3, sig_otu_stools_day3)
otu_day4_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 4, sig_otu_stools_day4)
otu_day5_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 5, sig_otu_stools_day5)
otu_day6_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 6, sig_otu_stools_day6)
otu_day7_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 7, sig_otu_stools_day7)
otu_day8_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 8, sig_otu_stools_day8)
otu_day9_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 9, sig_otu_stools_day9)
otu_day10_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 10, sig_otu_stools_day10)
otu_day15_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 15, sig_otu_stools_day15)
otu_day30_stats_stools <- pairwise_day_otu(all_otu_stools, "stools", 30, sig_otu_stools_day30)

#Think about how this code could be more DRY
#Combine pairwise dataframes for all days
otu_pairwise_stools <- rbind(otu_dayn5_stats_stools, otu_dayn1_stats_stools, otu_day0_stats_stools,
                             otu_day2_stats_stools, otu_day3_stats_stools, otu_day4_stats_stools,
                             otu_day5_stats_stools, otu_day6_stats_stools, otu_day7_stats_stools,
                             otu_day8_stats_stools, otu_day9_stats_stools, otu_day10_stats_stools,
                             otu_day15_stats_stools, otu_day30_stats_stools)

#Examine OTUs that are different between WM & M1 mice
WM_M1_otu_stools <- otu_pairwise_stools %>% 
  filter(group1 %in% c("WM", "M1") & group2 %in% c("WM", "M1")) %>% 
  filter(p.adj < 0.05)
WM_1RM1_otu_stools <- otu_pairwise_stools %>% 
  filter(group1 %in% c("WM", "1RM1") & group2 %in% c("WM", "1RM1")) %>% 
  filter(p.adj < 0.05)
#Lachnospiraceae and Porphyromonadaceae seem to be most of the OTUs that differ

#Examine OTUs that differ from either M1, 1RM1, or C groups) at day 5
d5_otu_stools <- otu_pairwise_stools %>% 
  filter(group1 %in% c("C", "M1", "1RM1") | group2 %in% c("C", "M1", "1RM1")) %>% 
  filter(p.adj < 0.05 & day == 5)

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

#Create heatmap of OTUs that differentiate PEG-3350 groups' stool samples-----
peg_sig_otus_1vs5 <- otu_pairwise_stools %>%
  filter(day >2) %>%  #Limit to when C. diff levels are starting to drop in 1-day PEG mice
  filter(group1 %in% c("WM", "M1", "1RM1") & group2 %in% c("WM", "M1", "1RM1")) %>% 
  filter(p.adj < 0.05) %>% 
  arrange(p.adj) %>% 
  distinct(otu) %>% 
  pull(otu)

peg_otus_M1_vs_FRM <- otu_pairwise_stools %>% 
  filter(day > 2) %>% #Limit to when C. diff levels are starting to drop in 1-day PEG mice
  filter(group1 %in% c("M1", "1RM1", "FRM") & group2 %in% c("M1", "1RM1", "FRM")) %>% 
  filter(p.adj < 0.05) %>% 
  arrange(p.adj) %>% 
  distinct(otu) %>% 
  pull(otu)

#PEG OTUs of interest: differentiate 1 vs 5-day PEG mice and 1-day vs Post_CDI PEG mice
peg_otus <- intersect_all(peg_sig_otus_1vs5, peg_otus_M1_vs_FRM)

#Heatmap of PEG OTUs of interest
peg_stools <- all_otu_stools %>% 
  filter(day %in% c("3", "4", "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Limit to when C. diff levels are starting to drop in 1-day PEG mice
  filter(group %in% c("WM", "WMC", "WMR",
                      "M1", "1RM1",
                      "CWM", "FRM", "RM","WMN")) 
peg_stool_days <- peg_stools %>% distinct(day) %>% pull(day)
facet_labels <- c( "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 + 10-day recovery",
                   "1-day PEG 3350", "1-day PEG 3350 + 1-day recovery",
                   "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350",
                   "5-day PEG 3350 without infection") #Create descriptive labels for facets
names(facet_labels) <- c("WM", "WMC", "WMR",
                         "M1", "1RM1",
                         "CWM", "FRM", "RM",
                         "WMN") #values that correspond to group, which is the variable we're faceting by
peg_stool <- hm_plot_otus(peg_stools, peg_otus[0:25], peg_stool_days)+
  scale_x_discrete(breaks = c(3:10, 15, 20, 25, 30), labels = c(3:10, 15, 20, 25, 30))
save_plot(filename = "results/figures/all_otus_heatmap_stools_peg.png", peg_stool, base_height = 15, base_width = 18)

#Heatmap of single OTUs of interest across all groups----
hm_1_otu(all_otu_stools, "Oscillibacter (OTU 45)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_oscillibacter_45.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Oscillibacter (OTU 57)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_oscillibacter_57.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Enterobacteriaceae (OTU 2)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_enterobacteriaceae_2.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Porphyromonadaceae (OTU 8)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_porphyromonadaceae_8.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lactobacillus (OTU 9)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lactobacillus_9.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Peptostreptococcaceae (OTU 12)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_peptostreptococcaceae_12.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 4)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_4.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 30)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_30.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 31)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_31.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 33)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_33.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 29)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_29.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 34)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_34.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 41)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_41.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Bacteria (OTU 43)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_bacteria_43.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Clostridium (OTU 51)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_clostridium_51.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Ruminococcaceae (OTU 50)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_ruminococcaceae_50.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Ruminococcaceae (OTU 52)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_ruminococcaceae_52.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Ruminococcaceae (OTU 54)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_ruminococcaceae_54.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Ruminococcaceae (OTU 60)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_ruminococcaceae_60.png", base_height = 15, base_width = 18)
ruminococcaceae_92 <- hm_1_otu(all_otu_stools, "Ruminococcaceae (OTU 92)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_ruminococcaceae_92.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Akkermansia (OTU 3)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_akkermansia_3.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Blautia (OTU 19)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_blautia_19.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Acetatifactor (OTU 55)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_acetatifactor_55.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Bacteria (OTU 108)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_bacteria_108.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Bacteroides (OTU 1)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_bacteroides_1.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Bifidobacterium (OTU 28)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_bifidobacterium_28.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Butyricicoccus (OTU 56)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_butyricicoccus_56.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Clostridium (OTU 22)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_clostridium_22.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 11)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_11.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 121)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_121.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 125)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_125.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 128)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_128.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 20)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_20.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 24)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_24.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 26)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_26.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 32)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_32.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 46)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_46.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 47)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_47.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 62)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_62.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 73)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_73.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lachnospiraceae (OTU 75)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lachnospiraceae_75.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Lactobacillus (OTU 13)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_lactobacillus_13.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Porphyromonadaceae (OTU 190)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_porphyromonadaceae_190.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Porphyromonadaceae (OTU 212)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_porphyromonadaceae_212.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Porphyromonadaceae (OTU 7)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_porphyromonadaceae_7.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Porphyromonadaceae (OTU 80)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_porphyromonadaceae_80.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Ruminococcaceae (OTU 59)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_ruminococcaceae_59.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Ruminococcaceae (OTU 98)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_ruminococcaceae_98.png", base_height = 15, base_width = 18)
hm_1_otu(all_otu_stools, "Turicibacter (OTU 5)", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_otus_heatmap_stools_turicibacter_5.png", base_height = 15, base_width = 18)

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

#Test for differences in genera relative abundances across groups at specific timepoints----
#subset agg_otu_data
all_genus <- agg_genus_data %>% 
  filter(group %in% all_groups)
all_genus_stools <- subset_stool(all_genus)
all_genus_tissues <- subset_tissue(all_genus)
#Function to test stool or tissue samples at the genus level:
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
  #Adjust p-values for testing multiple genera
  genus_stats_adjust <- genus_stats %>%
    select(-data) %>% #Keep everything but the data column
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
    arrange(p.value.adj) %>%
    write_tsv(path = paste0("data/process/all_genus_stats_day_", timepoint, "_", sample_type, ".tsv"))
}

#Create empty data frame to combine stat dataframes for all days that were tested
kw_genus_stools <- data.frame(genus=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                            C =double(), WM =double(),WMC=double(),WMR =double(),M1 =double(),`1RM1` =double(),
                            CWM =double(),FRM=double(),RM =double(),
                            CN=double(),WMN =double(),
                            p.value.adj=double(),day=double())%>% 
  rename(`1RM1` = X1RM1) #Rename 1RM1 to get rid of X

# Perform kruskal wallis tests at the genus level for all mice stool samples----
for (d in stool_test_days){
  kruskal_wallis_genus(d, all_genus_stools, "stools")
  #Make a list of significant genera across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/all_genus_stats_day_", d, "_stools.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_genus_stools_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, genus))
  kw_genus_stools <- add_row(kw_genus_stools, stats)  #combine all the dataframes together
}

# Perform pairwise Wilcoxan rank sum tests for genera that were significantly different across sources of mice on a series of days----
pairwise_day_genus <- function(sample_df, sample_type, timepoint, sig_genus_dayX){
  genus_stats <- sample_df %>% 
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
    write_tsv(path = paste0("data/process/all_genus_stats_day_", timepoint, "_", sample_type, "_sig.tsv"))
  #Format pairwise stats to use with ggpubr package
  plot_format_stats <- pairwise_stats %>% 
    select(-method) %>% 
    group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
    lapply(tidy_pairwise_genus) %>% 
    bind_rows() %>% 
    filter(!is.na(group2)) %>% #Remove rows that don't have a 2nd comparison group
    arrange(p.adj) %>% #Arrange by adjusted p value column 
    mutate(day = timepoint)
  return(plot_format_stats)  
}

#Pairwise Tests of genera that were significant on each day tested
#List of significant days:
pairwise_days_stools <- kw_genus_stools %>% filter(p.value.adj < 0.05) %>% 
  distinct(day) %>% pull(day)
#Pairwise test of significant days and their corresponding genera for stool samples 
#Create empty placeholder data frame
pairwise_genus_stools <- data.frame(genus=character(), group1=character(), group2=character(),
                                  p.adj = double(), day= double())
genus_dayn5_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", -5, `sig_genus_stools_day-5`)
genus_dayn1_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", -1, `sig_genus_stools_day-1`)
genus_day0_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 0, sig_genus_stools_day0)
genus_day1_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 1, sig_genus_stools_day1)
genus_day2_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 2, sig_genus_stools_day2)
genus_day3_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 3, sig_genus_stools_day3)
genus_day4_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 4, sig_genus_stools_day4)
genus_day5_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 5, sig_genus_stools_day5)
genus_day6_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 6, sig_genus_stools_day6)
genus_day7_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 7, sig_genus_stools_day7)
genus_day8_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 8, sig_genus_stools_day8)
genus_day9_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 9, sig_genus_stools_day9)
genus_day10_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 10, sig_genus_stools_day10)
genus_day15_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 15, sig_genus_stools_day15)
genus_day30_stats_stools <- pairwise_day_genus(all_genus_stools, "stools", 30, sig_genus_stools_day30)

#Think about how this code could be more DRY
#Combine pairwise dataframes for all days
genus_pairwise_stools <- rbind(genus_dayn5_stats_stools, genus_dayn1_stats_stools, genus_day0_stats_stools,
                             genus_day2_stats_stools, genus_day3_stats_stools, genus_day4_stats_stools,
                             genus_day5_stats_stools, genus_day6_stats_stools, genus_day7_stats_stools,
                             genus_day8_stats_stools, genus_day9_stats_stools, genus_day10_stats_stools,
                             genus_day15_stats_stools, genus_day30_stats_stools)

#Examine genera that are different between WM & M1 mice
WM_M1_genus_stools <- genus_pairwise_stools %>% 
  filter(group1 %in% c("WM", "M1") & group2 %in% c("WM", "M1")) %>% 
  filter(p.adj < 0.05)
WM_1RM1_genus_stools <- genus_pairwise_stools %>% 
  filter(group1 %in% c("WM", "1RM1") & group2 %in% c("WM", "1RM1")) %>% 
  filter(p.adj < 0.05)

#Examine genera that are different between RM & FRM mice (These likely are unrelated to C. difficile clearance)
RM_FRM_genus_stools <- genus_pairwise_stools %>% 
  filter(group1 %in% c("RM", "FRM") & group2 %in% c("RM", "FRM")) %>% 
  filter(p.adj < 0.05)
#Porpyhromonadaceae across 6 timepoints & Bacteroidales across 5 timepoints

#Examine genera that are different on day 5 between either C, M1, 1RM1 mice that clear within 10 days and WM, WMC, WMR, CWM, RM, FRM mice that remained colonized or clear on day 15
clear_v_prolonged_genus_stools <- genus_pairwise_stools %>% 
  filter(day == 5) %>% 
  filter(group1 %in% c("C", "M1", "1RM1") | group1 %in% c("WM", "WMC", "WMR", "CWM", "RM", "FRM"),
         group2 %in% c("C", "M1", "1RM1") | group2 %in% c("WM", "WMC", "WMR", "CWM", "RM", "FRM")) %>% 
  filter(p.adj < 0.05) %>% 
  group_by(genus) %>% 
  count() %>% 
  filter(n >= 15)

#Examine genera that are different on day 5 between PEG-treated mice: either M1, 1RM1 mice that clear within 10 days and WM, WMC, WMR, CWM, RM, FRM mice that remained colonized or clear on day 15
peg_clear_v_prolonged_genus_stools <- genus_pairwise_stools %>% 
  filter(day == 5) %>% 
  filter(group1 %in% c("M1", "1RM1") | group1 %in% c("WM", "WMC", "WMR", "CWM", "RM", "FRM"),
         group2 %in% c("M1", "1RM1") | group2 %in% c("WM", "WMC", "WMR", "CWM", "RM", "FRM")) %>% 
  filter(p.adj < 0.05) %>% 
  group_by(genus) %>% 
  count() %>% 
  pull(genus)

#Make an area plot of the genera that vary between PEG mice that clear within 10 days and PEG mice that remain colonized----
facet_labels <- all_labels <- c( "Clind.", 
                                 "1-day PEG 3350 + 1-day recovery", "1-day PEG 3350", 
                                 "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", 
                                 "Clind. + 3-day recovery + 1-day PEG 3350", "Clind. + 1-day PEG 3350",
                                 "5-day PEG 3350 + 10-day recovery",
                                 "5-day PEG 3350 + Clind.", 
                                 "5-day PEG 3350"
                                 )
#Create descriptive labels for facets
names(facet_labels) <- c("C", "1RM1", "M1", "FRM", "RM", "CWM", "WMR", "WMC", "WM") #values that correspond to group, which is the variable we're faceting by
area_plot <- all_genus_stools %>%
  #Solve different baseline timepoints across groups by specifying baseline for each group
  mutate(day = case_when(group == "C" & day %in% c("-15", "-11", "-1") ~ "B",
                         group == "WM" & day == "-5" ~ "B",
                         group == "WMC" & day == "-5" ~ "B",
                         group == "WMR" & day == "-15" ~ "B",
                         group == "CWM" & day == "-1" ~ "B",
                         group == "RM" & day == "-1" ~ "B",
                         group == "FRM" & day == "-1" ~ "B",
                         group == "M1" & day == "-11" ~ "B",
                         group == "1RM1" & day == "-2" ~ "B",
                         TRUE ~ day)) %>% 
  filter(day %in% c("B", "0", "1", "2", "3", "4", "5", "6", "7", "10")) %>% 
  filter(group %in% c("C", "1RM1", "M1", "FRM", "RM", "CWM", "WMR", "WMC", "WM")) %>% 
  mutate(group = fct_relevel(group, "C", "1RM1", "M1", "FRM", "RM", "CWM", "WMR", "WMC","WM")) %>% #Specify the order of the groups
  mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
  mutate(day = fct_relevel(day, "B", "0", "1", "2", "3", "4",
                           "5", "6", "7", "10")) %>% 
  filter(genus %in% peg_clear_v_prolonged_genus_stools) %>%
  group_by(group, day, genus) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  ggplot()+
  geom_area(aes(x = day, y=median, group = genus, fill = genus))+
  scale_fill_viridis(discrete = T)+
  labs(title=NULL,
       x=NULL,
       y=NULL)+
  facet_wrap(~group, labeller = labeller(group = facet_labels))+
  theme_classic()+
  theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
        legend.text = element_text(face = "italic"), #Italicize genus name
        legend.title = element_blank(),
        text = element_text(size = 16)) # Change font size for entire plot
save_plot(filename = "results/figures/all_genus_area_stools.png", area_plot, base_height = 15, base_width = 18)
ggplotly(area_plot)

#Examine the genera that vary between PEG mice that clear within 10 days and PEG mice that remain colonized in the tissue samples----
facet_labels <- all_labels <- c("cecum", "proximal colon", "distal colon",
                                "Clind.", "5-day PEG 3350",
                                  "5-day PEG 3350 + 10-day recovery")
#Create descriptive labels for facets
names(facet_labels) <- c("cecum", "proximal_colon", "distal_colon", 
                         "C", "WM", "WMR") #values that correspond to group, which is the variable we're faceting by
area_plot_tissues <- all_genus_tissues  %>%
  filter(day %in% c("4", "6", "20", "30")) %>% 
  filter(group %in% c("C", "WM", "WMR")) %>% 
  mutate(group = fct_relevel(group, "C", "WM", "WMR")) %>% #Specify the order of the groups
  mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
  mutate(day = fct_relevel(day, "4", "6", "20", "30")) %>% 
  mutate(sample_type = fct_relevel(sample_type, "cecum", "proximal_colon", "distal_colon")) %>% 
  filter(genus %in% peg_clear_v_prolonged_genus_stools) %>%
  group_by(group, day, genus, sample_type) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  ggplot()+
  geom_area(aes(x = day, y=median, group = genus, fill = genus))+
  scale_fill_viridis(discrete = T)+
  labs(title=NULL,
       x=NULL,
       y=NULL)+
  facet_grid(group ~ sample_type, labeller = labeller(group = facet_labels, sample_type = facet_labels))+
  theme_classic()+
  theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
        legend.text = element_text(face = "italic"), #Italicize genus name
        legend.title = element_blank(),
        text = element_text(size = 16)) # Change font size for entire plot
save_plot(filename = "results/figures/all_genus_area_tissues.png", area_plot_tissues, base_height = 15, base_width = 18)
ggplotly(area_plot_tissues)

#Examine genera that differ from either M1, 1RM1, or C groups) at day 5
d5_genus_stools <- genus_pairwise_stools %>% 
  filter(group1 %in% c("C", "M1", "1RM1") | group2 %in% c("C", "M1", "1RM1")) %>% 
  filter(p.adj < 0.05 & day == 5) %>% 


#Create empty data frame to combine stat dataframes for all days that were tested
kw_genus_tissues <- data.frame(genus=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                             C =double(), WM =double(),WMC=double(),WMR =double(),
                             CWM =double(),
                             CN=double(),WMN =double(),
                             p.value.adj=double(),day=double())
#Create empty data frames to combine stat dataframes for all days that were tested
kw_genus_cecum <- data.frame(genus=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                               C =double(), WM =double(),WMC=double(),WMR =double(),
                               CWM =double(),
                               CN=double(),WMN =double(),
                               p.value.adj=double(),day=double())
kw_genus_pc <- data.frame(genus=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                             C =double(), WM =double(),WMC=double(),WMR =double(),
                             CWM =double(),
                             CN=double(),WMN =double(),
                             p.value.adj=double(),day=double())
kw_genus_dc <- data.frame(genus=character(), statistic=double(), p.value = double(), parameter=double(), method=character(),
                             C =double(), WM =double(),WMC=double(),WMR =double(),
                             CWM =double(),
                             CN=double(),WMN =double(),
                             p.value.adj=double(),day=double())

# Perform kruskal wallis tests at the genus level for the tissue samples----
#do separate tests for each type of tissue
all_genus_cecum <- all_genus_tissues %>% filter(sample_type == "cecum")
all_genus_pc <- all_genus_tissues %>% filter(sample_type == "proximal_colon")
all_genus_dc <- all_genus_tissues %>% filter(sample_type == "distal_colon")
#for (d in tissue_test_days){
#  kruskal_wallis_genus(d, all_genus_tissues, "tissues")
#  #Make a list of significant genera across sources of mice for a specific day
#  stats <- read_tsv(file = paste0("data/process/all_genus_stats_day_", d, "_tissues.tsv")) %>% 
#    mutate(day = d)#Add a day column to specify day tested
#  name <- paste("sig_genus_tissues_day", d, sep = "")
#  assign(name, pull_significant_taxa(stats, genus))
#  kw_genus_tissues <- add_row(kw_genus_tissues, stats)  #combine all the dataframes together
}

for (d in tissue_test_days){
  kruskal_wallis_genus(d, all_genus_cecum, "cecum")
  #Make a list of significant genera across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/all_genus_stats_day_", d, "_cecum.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_genus_cecum_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, genus))
  kw_genus_cecum <- add_row(kw_genus_cecum, stats)  #combine all the dataframes together
}
for (d in tissue_test_days){
  kruskal_wallis_genus(d, all_genus_pc, "pc")
  #Make a list of significant genera across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/all_genus_stats_day_", d, "_pc.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_genus_pc_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, genus))
  kw_genus_pc <- add_row(kw_genus_pc, stats)  #combine all the dataframes together
}
for (d in tissue_test_days){
  kruskal_wallis_genus(d, all_genus_dc, "dc")
  #Make a list of significant genera across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/all_genus_stats_day_", d, "_dc.tsv")) %>% 
    mutate(day = d)#Add a day column to specify day tested
  name <- paste("sig_genus_dc_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, genus))
  kw_genus_dc <- add_row(kw_genus_dc, stats)  #combine all the dataframes together
}
kw_genus_cecum <- kw_genus_cecum %>% mutate(sample_type = "cecum")
kw_genus_pc <- kw_genus_pc %>% mutate(sample_type = "proximal_colon")
kw_genus_dc <- kw_genus_dc %>% mutate(sample_type = "distal_colon")
kw_genus_tissues <- kw_genus_cecum %>% add_row(kw_genus_pc) %>% add_row(kw_genus_dc)
#Pairwise Tests of genera that were significant on each day tested
#List of significant days:
pairwise_days_tissues <- kw_genus_tissues %>% filter(p.value.adj < 0.05) %>% 
  distinct(day) %>% pull(day)
pairwise_days_cecum <- kw_genus_cecum %>% filter(p.value.adj < 0.05) %>% 
  distinct(day) %>% pull(day)
#Pairwise test of significant days and their corresponding genera for stool samples 
#Create empty placeholder data frame
pairwise_genus_tissues <- data.frame(genus=character(), group1=character(), group2=character(),
                                    p.adj = double(), day= double())
pairwise_genus_cecum <- data.frame(genus=character(), group1=character(), group2=character(),
                                     p.adj = double(), day= double())
genus_day4_stats_cecum <- pairwise_day_genus(all_genus_cecum, "cecum", 4, sig_genus_cecum_day4)
genus_day6_stats_cecum <- pairwise_day_genus(all_genus_cecum, "cecum", 6, sig_genus_cecum_day6)
genus_day30_stats_cecum <- pairwise_day_genus(all_genus_cecum, "cecum", 30, sig_genus_cecum_day30)

#Overlapping genera over time in the cecum:
overlap_cecum <- intersect_all(sig_genus_cecum_day4, sig_genus_cecum_day6, sig_genus_cecum_day30)

#Think about how this code could be more DRY
#Combine pairwise dataframes for all days
genus_pairwise_stools <- rbind(genus_dayn5_stats_stools, genus_dayn1_stats_stools, genus_day0_stats_stools,
                               genus_day2_stats_stools, genus_day3_stats_stools, genus_day4_stats_stools,
                               genus_day5_stats_stools, genus_day6_stats_stools, genus_day7_stats_stools,
                               genus_day8_stats_stools, genus_day9_stats_stools, genus_day10_stats_stools,
                               genus_day15_stats_stools, genus_day30_stats_stools)


#Create heatmap of significant genera for all stool samples----
#Rank genera by adjusted p-value
hm_sig_genus_p_adj <- kw_genus_stools %>% 
  filter(p.value.adj < 0.05) %>% 
  arrange(p.value.adj) %>% 
  distinct(genus) %>% 
  slice_head(n = 25) %>% 
  pull(genus)

hm_stool_days <- all_diversity_stools %>% distinct(day) %>% pull(day)
facet_labels <- all_labels #Create descriptive labels for facets
names(facet_labels) <- c("C", "WM", "WMC", "WMR",
                         "M1", "1RM1",
                         "CWM", "FRM", "RM",
                         "CN", "WMN") #values that correspond to group, which is the variable we're faceting by
hm_stool <- hm_plot_genus(all_genus_stools, hm_sig_genus_p_adj, hm_stool_days)+
  scale_x_discrete(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30), labels = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30)) 
save_plot(filename = "results/figures/all_genus_heatmap_stools.png", hm_stool, base_height = 15, base_width = 18)

#Heatmap of single OTUs of interest across all groups----
hm_1_genus(all_genus_stools, "Phenylobacterium", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_phenylobacterium.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Enterobacteriaceae Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_enterobacteriaceae.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Bacteroides", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_bacteroides.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Peptostreptococcaceae Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_peptostreptococcaceae.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Ruminococcaceae Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_ruminococcaceae.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Bifidobacterium", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_bifidobacterium.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Lactobacillus", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_lactobacillus.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Porphyromonadaceae Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_porphyromonadaceae.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Clostridiales Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_clostridiales.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Oscillibacter", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_oscillibacter.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Bacteroidales Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_bacteroidales.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_unclassified.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Blautia", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_blautia.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Lachnospiraceae Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_lachnospiraceae.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Acetatifactor", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_acetatifactor.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Butyricicoccus", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_butyricicoccus.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Clostridium XlVb", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_clostridium_XlVb.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Turicibacter", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_turicibacter.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Firmicutes Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_firmicutes.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Olsenella", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_olsenella.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Alistipes", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_alistipes.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Chitinophagaceae Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_chitinophagaceae.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Erysipelotrichaceae Unclassified", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_erysipelotrichaceae.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Akkermansia", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_akkermansia.png", base_height = 15, base_width = 18)
hm_1_genus(all_genus_stools, "Clostridium IV", hm_stool_days) %>% 
  save_plot(filename = "results/figures/all_genus_heatmap_stools_clostridium_IV.png", base_height = 15, base_width = 18)

#Examine correlation between C. difficile CFUs and Peptostreptococcaceae (OTU 12) relative abundance-----
pepto_cfu_otu <- all_otu_stools %>% 
  filter(otu == "Peptostreptococcaceae (OTU 12)") %>% 
  filter(!is.na(avg_cfu)) 

pepto_cfu_otu_plot <- pepto_cfu_otu%>% 
  ggplot(aes(x = avg_cfu, y = agg_rel_abund))+
  geom_point()+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  scale_x_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12))+ #Scientific notation labels for y-axis
  theme_classic()+
  labs(x = "CFU/g Feces",
       y = "Relative abundance (%)")+
  geom_smooth(method = "lm")
ggsave("exploratory/notebook/c_diff_cfu_otu_correlation.png", pepto_cfu_otu_plot)

#Examine correlation 
cor.test(pepto_cfu_otu$avg_cfu, pepto_cfu_otu$agg_rel_abund, method = "spearman")
#rho = 0.79, p-value 2.2e-16   

#Examine peptostreptococcoaceae over time in the stools of all mice with 16S sequencing data
pepto_cfu <- all_otu_stools %>% 
  filter(otu == "Peptostreptococcaceae (OTU 12)") %>% 
  mutate(day = as.integer(day)) #Fix day variable
pepto_otu_name <- pepto_cfu %>% 
  pull(otu_name)
pepto_otu_median <- pepto_cfu %>% 
  group_by(group, day) %>% 
  summarize(median=(median(agg_rel_abund + 1/2000))) %>% 
  ungroup
pepto_otu_mice <- pepto_cfu %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>%
  select(day, agg_rel_abund, otu, group)
pepto_over_time <-  ggplot(NULL)+
  geom_point(pepto_otu_mice, mapping = aes(x=day, y=agg_rel_abund), size  = 1.5, position = position_dodge(width = 0.6))+
  geom_line(pepto_otu_median, mapping = aes(x=day, y=median, group = group), size = 1, show.legend = FALSE)+
  geom_hline(yintercept=1/1000, color="gray")+
  labs(title=pepto_otu_name,
       x="Day",
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme_classic()+
  theme(plot.title=element_markdown(hjust = 0.5),
        panel.grid.minor.x = element_line(size = 0.4, color = "grey"),  # Add gray lines to clearly separate symbols by days)
        text = element_text(size = 18)) # Change font size for entire plot
ggsave("exploratory/notebook/c_diff_otu_time_all_mice.png", pepto_over_time, height = 4, width = 8.5)

#Examine Peptostreptococcaceae OTUs in the tissues----
pepto_cfu_otu_tissues <- all_otu_tissues %>% 
  filter(otu == "Peptostreptococcaceae (OTU 12)") %>% 
  filter(!is.na(avg_cfu)) 

pepto_cfu_otu_plot_tissues <- pepto_cfu_otu_tissues %>% 
  ggplot(aes(x = avg_cfu, y = agg_rel_abund))+
  geom_point()+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  scale_x_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12))+ #Scientific notation labels for y-axis
  theme_classic()+
  labs(x = "CFU/g Feces",
       y = "Relative abundance (%)")+
  geom_smooth(method = "lm")
ggsave("exploratory/notebook/c_diff_cfu_otu_correlation_tissues.png", pepto_cfu_otu_plot_tissues)

#Examine correlation 
cor.test(pepto_cfu_otu_tissues$avg_cfu, pepto_cfu_otu_tissues$agg_rel_abund, method = "spearman")
#rho = 0.77, p-value 2.2e-16  

#Examine correlations across tissue types and stools
pepto_cfu_otu_all <- all_otu %>% 
  filter(sample_type %in% c("cecum", "proximal_colon", "distal_colon", "stool")) %>% 
  filter(otu == "Peptostreptococcaceae (OTU 12)") %>% 
  filter(!is.na(avg_cfu))
  
facet_labels <- c("Cecum", "Proximal colon", "Distal Colon", "Stool") #Create descriptive labels for facets
names(facet_labels) <- c("cecum", "proximal_colon", "distal_colon", "stool") #values that correspond to group, which is the variable we're faceting by
pepto_cfu_otu_plot_sample_type <- pepto_cfu_otu_all %>% 
  mutate(sample_type = fct_relevel(sample_type, "cecum", "proximal_colon", "distal_colon", "stool")) %>% 
  ggplot(aes(x = avg_cfu, y = agg_rel_abund))+
  geom_point()+
  facet_wrap(~sample_type, labeller = labeller(sample_type = facet_labels))+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  scale_x_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12))+ #Scientific notation labels for y-axis
  theme_classic()+
  labs(x = "CFU/g Feces",
       y = "Relative abundance (%)")+
  geom_smooth(method = "lm")+
  theme(strip.background = element_blank()) #get rid of box around facet_wrap labels
ggsave("exploratory/notebook/c_diff_cfu_otu_correlation_sample_type.png", pepto_cfu_otu_plot_sample_type)

#Examine correlations across tissue types----
pepto_cecum <- pepto_cfu_otu_tissues %>% filter(sample_type == "cecum")
pepto_pc <- pepto_cfu_otu_tissues %>% filter(sample_type == "proximal_colon")
pepto_dc <- pepto_cfu_otu_tissues %>% filter(sample_type == "distal_colon")

#Examine correlation between plated C. diff, and Peptostreptococcace unclassified in the cecum
cor.test(pepto_cecum$avg_cfu, pepto_cecum$agg_rel_abund, method = "spearman")
#rho = 0.85, p-value 2.2e-16  
#Examine correlation between plated C. diff, and Peptostreptococcace unclassified in the proximal colon
cor.test(pepto_pc$avg_cfu, pepto_pc$agg_rel_abund, method = "spearman")
#rho = 0.80, p-value 2.2e-16 
#Examine correlation between plated C. diff, and Peptostreptococcace unclassified in the distal colon
cor.test(pepto_dc$avg_cfu, pepto_dc$agg_rel_abund, method = "spearman")
#rho = 0.64, p-value 2.2e-16 

#Examine peptostreptococcoaceae over time in the stools of all mice with 16S sequencing data
pepto_cfu_t <- all_otu_tissues %>% 
  filter(otu == "Peptostreptococcaceae (OTU 12)") %>% 
  mutate(day = as.integer(day)) #Fix day variable
pepto_otu_name_t <- pepto_cfu_t %>% 
  pull(otu_name)
pepto_otu_median_t <- pepto_cfu_t %>% 
  group_by(group, day) %>% 
  summarize(median=(median(agg_rel_abund + 1/2000))) %>% 
  ungroup
pepto_otu_mice_t <- pepto_cfu_t %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>%
  select(day, agg_rel_abund, otu, group)
pepto_over_time_t <-  ggplot(NULL)+
  geom_point(pepto_otu_mice_t, mapping = aes(x=day, y=agg_rel_abund), size  = 1.5, position = position_dodge(width = 0.6))+
  geom_line(pepto_otu_median_t, mapping = aes(x=day, y=median, group = group), size = 1, show.legend = FALSE)+
  geom_hline(yintercept=1/1000, color="gray")+
  labs(title=pepto_otu_name,
       x="Day",
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  scale_x_continuous(breaks = c(-15, -10, -5, -4, -2, -1:10, 15, 20, 30),
                     limits = c(-16,31),
                     minor_breaks = c(-15.5,-14.5, -10.5, -9.5, -5.5, -4.5, -3.5, -2.5, -1.5:10.5, 14.5, 15.5, 19.5, 20.5, 29.5, 30.5)) +
  theme_classic()+
  theme(plot.title=element_markdown(hjust = 0.5),
        panel.grid.minor.x = element_line(size = 0.4, color = "grey"),  # Add gray lines to clearly separate symbols by days)
        text = element_text(size = 18)) # Change font size for entire plot
ggsave("exploratory/notebook/c_diff_otu_time_all_mice_tissues.png", pepto_over_time_t, height = 4, width = 8.5)
