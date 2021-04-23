source("code/utilities.R") #Loads libraries, reads in metadata, functions
library(viridis)
library(scales)

set.seed(19760620) #Same seed used for mothur analysis

#Used to figure out color scale:
values = brewer_pal("qual", palette = "Dark2")(3)

#Function to read in model performances:
#file_path = path to performance results file
#taxa_level_name = "genus" or "otu"
#comp_name = name of comparison in quotes
read_perf_results <- function(file_path, taxa_level_name){
  read_csv(file_path) %>% 
    mutate(taxa_level = taxa_level_name) 
}
#Read in genus level model performance:
genus_results <- read_perf_results("results/performance_results.csv", 
                              "genus") 

#Examine AUC results for each machine learning method----
auc_results <- genus_results %>% 
  group_by(method) %>% 
  summarize(median = median(AUC)) %>% 
  arrange(desc(median))

#Plot performance for all methods with genus input data----
performance_genus <- genus_results %>% 
  filter(taxa_level == "genus") %>% 
  mutate(method = fct_relevel(method, c("rf", "glmnet", "svmRadial"))) %>% #Reorder methods so left to right is in order of descending AUC
  ggplot(aes(x = method, y = AUC, color = method)) +
  geom_boxplot(alpha=0.5, fatten = 4) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  scale_color_manual(values = c("#D95F02", "#1B9E77", "#E7298A"),                     
                     breaks=c("rf", "glmnet", "svmRadial"), 
                     labels = c("random forest", "logistic regression", "support vector machine")) +
  scale_y_continuous(name = "AUC",
                     breaks = seq(0.4, 1, 0.1),
                     limits=c(0.4, 1),
                     expand=c(0,0)) +
  labs(x = NULL)+
  scale_x_discrete(label = c("random forest", "logistic regression", "support vector machine"))+
  theme_bw()  +
  theme(legend.position = "none",
        text = element_text(size = 19),# Change font size for entire plot
#        axis.ticks.x = element_blank(), #Remove x axis ticks
#        axis.text.x = element_blank(), #Remove x axis text (all models are at the genus level),
        axis.text.x = element_text(angle = 45, hjust = 1), #Angle axis labels
        strip.background = element_blank()) +#Make Strip backgrounds blank
  guides(color=guide_legend(nrow = 2))+   #Legend in 2 rows so it doesn't get cut off 
  ggsave("results/figures/ml_performance_genus.png", height = 5, width = 8)

#Examine permutation importance results for random forest---- 
# random forest  had the highest median AUC
#Function to read in feature importance from random forest model----
#Will also get the corresponding genus name from taxonomy file
#And format genus name to work with glue package for italicizing labels
#file_path = path to the file name
read_feat_imp <- function(file_path){
  feat_imp <- read_csv(file_path)
  taxa_info <- read.delim('data/process/peg3350.genus.taxonomy', header=T, sep='\t') %>%
    select(-Size) %>%
    mutate(names=OTU) %>%
    select(-OTU)
  final_feat_imp <- inner_join(feat_imp, taxa_info, by="names") %>%
    ungroup() %>%    
    mutate(names=str_to_upper(names)) %>%
    mutate(genus=gsub("(.*);.*","\\1",Taxonomy)) %>%
    mutate(genus=gsub("(.*)_.*","\\1",Taxonomy)) %>%
    mutate(genus=gsub("(.*);.*","\\1",Taxonomy)) %>% 
    mutate(genus=str_replace_all(genus, c("_" = " ", #Removes all other underscores
                                          "unclassified" = "Unclassified"))) %>% 
    select(-Taxonomy)
  return(final_feat_imp)
}

rf_feat <- read_feat_imp("results/combined_feature-importance_rf.csv") #%>% 
  #Transform bactname variable into factor 
  #mutate(bactname = factor(bactname, levels = unique(as.factor(bactname))))

#Function to get the top 20 features that have the largest impact on AUROC
#df = dataframe of feature importances for the 100 seeds
top_20 <- function(df){
  data_first_20 <- df %>% 
    group_by(genus) %>% 
    summarize(median = median(perf_metric_diff)) %>% #Get the median performance metric diff. for each feature
    arrange(desc(median)) %>% #Arrange from largest median to smallest
    head(20)
  
  return(data_first_20)
}

#Top 20 genera for each input dataset with the random forest model
rf_top_feat <- top_20(rf_feat) %>% pull(genus)

#Function to filter to top genera for each pairwise comparison & plot results----
#df = dataframes of feature importances for all seeds
#top_feat = dataframes of top features
#comp_name = name of comparison to title the plot (in quotes)
plot_feat_imp <- function(df, top_feat){
  df %>% 
    filter(genus %in% top_feat) %>% 
    ggplot(aes(fct_reorder(genus, -perf_metric_diff, .desc = TRUE), perf_metric_diff, color = genus))+
    stat_summary(fun = 'median', 
                 fun.max = function(x) quantile(x, 0.75), 
                 fun.min = function(x) quantile(x, 0.25),
                 position = position_dodge(width = 1)) + 
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_bact,
                        labels=legend_bact)+
    coord_flip()+
    labs(title=NULL, 
         x=NULL,
         y="Difference in AUROC")+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          text = element_text(size = 15),# Change font size for entire plot
          axis.text.y = element_text(face = "italic"), #Genus in italics
          strip.background = element_blank(),
          legend.position = "none")   
}

#Plot feature importances for the top genera for each comparison----
#Figure out colors for top 5 most abundant genera based on viridis scale
show_col(viridis_pal()(5))
color_scheme_df <- rf_feat %>% 
  distinct(genus) %>% #Limit to unique genera
  filter(genus %in% rf_top_feat) %>% 
  #Assign viridis colors to top 5 most abundant genera
  mutate(color = case_when(#genus == "Porphyromonadaceae Unclassified" ~ "#440154FF",
                           #genus == "Akkermansia" ~ "#3B528BFF",
                           #genus == "Enterobacteriaceae Unclassified" ~ "#21908CFF",
                           #genus == "Lachnospiraceae Unclassified" ~ "#5DC863FF",
                           #genus == "Bacteroides" ~ "#FDE725FF",
                           TRUE ~ "black")) #Rest of colors should be black
color_scheme <- color_scheme_df %>% pull(color)
color_bact <- color_scheme_df %>% pull(genus)
legend_bact <- color_scheme_df %>% pull(genus)
rf_feat_5dpi <- plot_feat_imp(rf_feat, rf_top_feat)+
  ggsave("results/figures/ml_top_features_genus.png", height = 5, width = 8)

#Examine relative abundances in mice that clear within 10 days vs mice with prolonged colonization----
source("code/16S_common_files.R") #Reads in mothur output files

#Create shape scale based on each subset group
shape_scheme <- c(1, 4, 19, 8)
shape_groups <- c("clind.", "1-day", "5-day", "post-CDI")
shape_labels <- c("Clind.", "1-day PEG", "5-day PEG", "Post-CDI PEG")

interp_genera_d5_top_10 <- rf_top_feat[1:10]
#Plot the top 10 features that were important to Day 5 model and 
#using facet_wrap, highlight the genera that correlate with colonization
top10_d5_model_taxa <- agg_genus_data %>% 
  filter(day == 5) %>% #Used d5 timepoint for ml input data
  filter(clearance_status_d10 %in% c("colonized", "cleared")) %>% #Remove samples we don't have clearance status d10 data for
  filter(genus %in% interp_genera_d5_top_10) %>%
  #Reorder genera to match contribution to model
  mutate(genus = fct_relevel(genus, interp_genera_d5_top_10)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% 
  group_by(clearance_status_d10, genus) %>% 
  mutate(median=(median(agg_rel_abund + 1/2000))) %>% #create a column of median values for each group
  ungroup() %>% 
  ggplot(aes(x=clearance_status_d10, y =agg_rel_abund, colour= clearance_status_d10))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median, ymin = median), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(aes(shape = subset), size=2, show.legend = TRUE, alpha = .4) +
  scale_colour_manual(name=NULL,
                      values=c("blue", "red"),
                      breaks=c("cleared", "colonized"),
                      labels=c("cleared", "colonized"))+
  scale_shape_manual(name=NULL,
                     values=shape_scheme,
                     breaks=shape_groups,
                     labels=shape_labels)+
  geom_hline(yintercept=1/1000, color="gray")+
  facet_wrap(~ genus, nrow=2, labeller = label_wrap_gen(width = 10))+
  labs(title=NULL,
       x=NULL,
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()+
  theme(text = element_text(size = 14),
        strip.text = element_text(hjust = 0.5, size = 8.6, face = "italic"),
        axis.text.x = element_blank(),
        legend.position = "bottom")
save_plot(filename = paste0("results/figures/ml_top10_d5_genus.png"), top10_d5_model_taxa, base_height = 5, base_width = 8)  

#Create area plot of genera that vary between mice that clear within 10 dpi and mice that have prolonged colonization
#List of genera to include in area plot, select genera where the median for either cleared or colonized mice is greater than 1%
interp_genera_d5_top_10_abundant <- c("Porphyromonadaceae Unclassified", "Akkermansia", "Enterobacteriaceae Unclassified",
                                      "Lachnospiraceae Unclassified", "Bacteroides")
#Labels for facet
facet_labels <- c("cleared", "colonized")
names(facet_labels) <- c("cleared", "colonized")
#Create plot
topabund_5_model_taxa_area_plot <- agg_genus_data %>% 
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
  filter(day %in% c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>%  #Use baseline through day 15 timepoints
  mutate(day = fct_relevel(day, "B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>% 
  filter(clearance_status_d10 %in% c("colonized", "cleared")) %>% #Remove samples we don't have clearance status d10 data for
  filter(genus %in% interp_genera_d5_top_10_abundant) %>%
  #Reorder genera to match contribution to model
  mutate(genus = fct_relevel(genus, interp_genera_d5_top_10_abundant)) %>% 
  group_by(clearance_status_d10, genus, day) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  ggplot()+
  geom_area(aes(x = day, y=median, group = genus, fill = genus))+
  scale_fill_viridis(discrete = T)+
  labs(title=NULL,
       x=NULL,
       y="Relative abundance (%)") +
#To do: Figure out how to get scale to match other plots  
#  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
  facet_wrap(~clearance_status_d10, labeller = labeller(clearance_status_d10 = facet_labels), nrow = 1)+
  theme_classic()+
  theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
        legend.text = element_text(face = "italic"), #Italicize genus name
        legend.title = element_blank(),
        text = element_text(size = 16))+ # Change font size for entire plot
  guides(color = guide_legend(ncol = 2))

#Extract legend for area plot
legend_area_plot <- get_legend(topabund_5_model_taxa_area_plot) %>% as_ggplot()
save_plot("results/figures/ml_abund_5_genus_area_legend.png", legend_area_plot, base_height = 1.5, base_width = 4)
#Save area plot without legend
topabund_5_model_taxa_area_plot <- topabund_5_model_taxa_area_plot+
  theme(legend.position = "none") #Remove legend
save_plot(filename = paste0("results/figures/ml_abund_5_genus_area.png"), topabund_5_model_taxa_area_plot, base_height = 4, base_width = 7)  

#Line plots instead of area plot to convey how the most abundant genera that were top contributors to machine learning models were changing over time
topabund_5_model_taxa_line_plot <- agg_genus_data %>% 
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
  filter(day %in% c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>%  #Use baseline through day 15 timepoints
  mutate(day = fct_relevel(day, "B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>% 
  filter(clearance_status_d10 %in% c("colonized", "cleared")) %>% #Remove samples we don't have clearance status d10 data for
  filter(genus %in% interp_genera_d5_top_10_abundant) %>%
  #Reorder genera to match contribution to model
  mutate(genus = fct_relevel(genus, interp_genera_d5_top_10_abundant)) %>% 
  group_by(clearance_status_d10, genus, day) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop")  %>% #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  ggplot()+
  geom_line(aes(x = day, y=median, color=clearance_status_d10, group = clearance_status_d10), linetype = "solid")+
  scale_colour_manual(name=NULL,
                      values=c("blue", "red"),
                      breaks=c("cleared", "colonized"),
                      labels=c("cleared", "colonized"))+
  scale_x_discrete(limits = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"), 
                   labels = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"),
                   breaks = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"))+
  scale_y_continuous(trans = "log10", limits = c(1/10900, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_hline(yintercept=1/1000, color="gray")+ #Represents limit of detection
  labs(title=NULL,
       x="Days Post-Infection",
       y="Relative abundance (%)")+
  facet_wrap(~genus, nrow = 1, labeller = label_wrap_gen(width = 12))+
  theme_classic()+
  theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
        panel.spacing = unit(1, "lines"), #Increase spacing between facets
        plot.title = element_markdown(hjust = 0.5), #Have only the genera names show up as italics
        text = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(face = "italic", size = 12),
        legend.position = "None")
save_plot(filename = paste0("results/figures/ml_abund_5_genus_lineplot.png"), topabund_5_model_taxa_line_plot, base_height = 3, base_width = 10)  

#Examine Porphyromondaceae OTUs that are similar to Muribaculum intestinale----
#See code/blast_otus.R for how blast search was performed and results
muribac_otus <- c("Porphyromonadaceae (OTU 6)", "Porphyromonadaceae (OTU 7)", "Porphyromonadaceae (OTU 8)", "Porphyromonadaceae (OTU 21)")
muribac_otu_names <- c("*Porphyromonadaceae* (OTU 6)", "*Porphyromonadaceae* (OTU 7)", "*Porphyromonadaceae* (OTU 8)", "*Porphyromonadaceae* (OTU 21)")
#Line plots instead of area plot to convey how the most abundant genera that were top contributors to machine learning models were changing over time
muribac_otus_line_plot <- agg_otu_data %>% 
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
  filter(day %in% c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>%  #Use baseline through day 15 timepoints
  mutate(day = fct_relevel(day, "B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>% 
  filter(clearance_status_d10 %in% c("colonized", "cleared")) %>% #Remove samples we don't have clearance status d10 data for
  filter(otu %in% muribac_otus) %>%
  #Reorder genera to match contribution to model
  mutate(otu_name = fct_relevel(otu_name, muribac_otu_names)) %>% #Reorganize otu_names to match OTU number
  group_by(clearance_status_d10, otu_name, day) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop")  %>% #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  ggplot()+
  geom_line(aes(x = day, y=median, color=clearance_status_d10, group = clearance_status_d10), linetype = "solid", alpha = 0.6)+
  scale_colour_manual(name=NULL,
                      values=c("blue", "red"),
                      breaks=c("cleared", "colonized"),
                      labels=c("cleared", "colonized"))+
  scale_x_discrete(limits = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"), 
                   labels = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"),
                   breaks = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"))+
  scale_y_continuous(trans = "log10", limits = c(1/10900, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_hline(yintercept=1/1000, color="gray")+ #Represents limit of detection
  labs(title=NULL,
       x="Days Post-Infection",
       y="Relative abundance (%)")+
  facet_wrap(~otu_name, nrow = 1, labeller = label_wrap_gen(width = 10))+
  theme_classic()+
  theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
        panel.spacing = unit(1, "lines"), #Increase spacing between facets
        plot.title = element_markdown(hjust = 0.5), #Have only the genera names show up as italics
        text = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        strip.text = element_markdown(size = 10),
        legend.position = "None")
save_plot(filename = paste0("results/figures/ml_abund_otu_muribaculum_lineplot.png"), muribac_otus_line_plot, base_height = 3, base_width = 10)  

#Examine other Porphyromondaceae OTUs, found by looking in taxonomy file----
porphyr_otus <- c("Porphyromonadaceae (OTU 10)", "Porphyromonadaceae (OTU 14)", "Porphyromonadaceae (OTU 15)", "Porphyromonadaceae (OTU 25)")
porphyr_otu_names <- c("*Porphyromonadaceae* (OTU 10)", "*Porphyromonadaceae* (OTU 14)", "*Porphyromonadaceae* (OTU 15)", "*Porphyromonadaceae* (OTU 25)")
#Line plots instead of area plot to convey how the most abundant genera that were top contributors to machine learning models were changing over time
porphyr_otus_line_plot <- agg_otu_data %>% 
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
  filter(day %in% c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>%  #Use baseline through day 15 timepoints
  mutate(day = fct_relevel(day, "B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>% 
  filter(clearance_status_d10 %in% c("colonized", "cleared")) %>% #Remove samples we don't have clearance status d10 data for
  filter(otu %in% porphyr_otus) %>%
  #Reorder genera to match contribution to model
  mutate(otu_name = fct_relevel(otu_name, porphyr_otu_names)) %>% #Reorganize otu_names to match OTU number
  group_by(clearance_status_d10, otu_name, day) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop")  %>% #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  ggplot()+
  geom_line(aes(x = day, y=median, color=clearance_status_d10, group = clearance_status_d10), linetype = "solid", alpha = 0.6)+
  scale_colour_manual(name=NULL,
                      values=c("blue", "red"),
                      breaks=c("cleared", "colonized"),
                      labels=c("cleared", "colonized"))+
  scale_x_discrete(limits = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"), 
                   labels = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"),
                   breaks = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"))+
  scale_y_continuous(trans = "log10", limits = c(1/10900, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_hline(yintercept=1/1000, color="gray")+ #Represents limit of detection
  labs(title=NULL,
       x="Days Post-Infection",
       y="Relative abundance (%)")+
  facet_wrap(~otu_name, nrow = 1, labeller = label_wrap_gen(width = 10))+
  theme_classic()+
  theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
        panel.spacing = unit(1, "lines"), #Increase spacing between facets
        plot.title = element_markdown(hjust = 0.5), #Have only the genera names show up as italics
        text = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        strip.text = element_markdown(size = 10),
        legend.position = "None")
save_plot(filename = paste0("results/figures/ml_abund_otu_porphyromonadaceae_lineplot.png"), porphyr_otus_line_plot, base_height = 3, base_width = 10)  

#Examine Lachnospiraceae OTUs, found by looking in taxonomy file----
lachno_otus <- c("Lachnospiraceae (OTU 4)", "Lachnospiraceae (OTU 11)", "Lachnospiraceae (OTU 16)", "Lachnospiraceae (OTU 17)")
lachno_otu_names <- c("*Lachnospiraceae* (OTU 4)", "*Lachnospiraceae* (OTU 11)", "*Lachnospiraceae* (OTU 16)", "*Lachnospiraceae* (OTU 17)")
#Line plots instead of area plot to convey how the most abundant genera that were top contributors to machine learning models were changing over time
lachno_otus_line_plot <- agg_otu_data %>% 
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
  filter(day %in% c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>%  #Use baseline through day 15 timepoints
  mutate(day = fct_relevel(day, "B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15")) %>% 
  filter(clearance_status_d10 %in% c("colonized", "cleared")) %>% #Remove samples we don't have clearance status d10 data for
  filter(otu %in% lachno_otus) %>%
  #Reorder genera to match contribution to model
  mutate(otu_name = fct_relevel(otu_name, lachno_otu_names)) %>% #Reorganize otu_names to match OTU number
  group_by(clearance_status_d10, otu_name, day) %>% 
  summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop")  %>% #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
  ggplot()+
  geom_line(aes(x = day, y=median, color=clearance_status_d10, group = clearance_status_d10), linetype = "solid", alpha = 0.6)+
  scale_colour_manual(name=NULL,
                      values=c("blue", "red"),
                      breaks=c("cleared", "colonized"),
                      labels=c("cleared", "colonized"))+
  scale_x_discrete(limits = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"), 
                   labels = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"),
                   breaks = c("B", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15"))+
  scale_y_continuous(trans = "log10", limits = c(1/10900, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  geom_hline(yintercept=1/1000, color="gray")+ #Represents limit of detection
  labs(title=NULL,
       x="Days Post-Infection",
       y="Relative abundance (%)")+
  facet_wrap(~otu_name, nrow = 1, labeller = label_wrap_gen(width = 10))+
  theme_classic()+
  theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
        panel.spacing = unit(1, "lines"), #Increase spacing between facets
        plot.title = element_markdown(hjust = 0.5), #Have only the genera names show up as italics
        text = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        strip.text = element_markdown(size = 10),
        legend.position = "None")
save_plot(filename = paste0("results/figures/ml_abund_otu_lachnospiraceae_lineplot.png"), lachno_otus_line_plot, base_height = 3, base_width = 10)  

  