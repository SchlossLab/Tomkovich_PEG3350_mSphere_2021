source("code/utilities.R") #Loads libraries, reads in metadata, functions

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
  ggsave("results/figures/ml_performance_otu.png", height = 5, width = 8)

#Examine permutation importance results for random forest---- 
# random forest  had the highest median AUC
#Function to read in feature importance from random forest model----
#Will also get the corresponding otu name from taxonomy file
#And format otu name to work with glue package for italicizing labels
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


#Function to filter to top OTUs for each pairwise comparison & plot results----
#df = dataframes of feature importances for all seeds
#top_otus = dataframes of top otus
#comp_name = name of comparison to title the plot (in quotes)
plot_feat_imp <- function(df, top_feat){
  df %>% 
    filter(genus %in% top_feat) %>% 
    ggplot(aes(fct_reorder(genus, -perf_metric_diff, .desc = TRUE), perf_metric_diff))+
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

#Plot feature importances for the top OTUs for each comparison----
rf_feat_5dpi <- plot_feat_imp(rf_feat, rf_top_feat)+
  ggsave("results/figures/ml_top_features_genus.png", height = 5, width = 8)

#Examine relative abundances in mice that clear within 10 days vs mice with prolonged colonization----
source("code/16S_common_files.R") #Reads in mothur output files

#Create shape scale based on each subset group
shape_scheme <- c(1, 4, 19, 8)
shape_groups <- c("clind.", "1-day", "5-day", "post-CDI")
shape_labels <- c("Clind.", "1-day", "5-day", "Post-CDI")

interp_genera_d5_top_10 <- rf_top_feat[1:10]
#Plot the top 10 features that were important to Day 5 model and 
#using facet_wrap, highlight the OTUs that correlate with colonization
top10_d5_model_taxa <- agg_genus_data %>% 
  filter(day == 5) %>% #Used d5 timepoint for ml input data
  filter(clearance_status_d10 %in% c("colonized", "cleared")) %>% #Remove samples we don't have clearance status d10 data for
  filter(genus %in% interp_genera_d5_top_10) %>%
  #Reorder genera to match contribution to model
  mutate(genus = fct_relevel(genus, interp_genera_d5_top_10)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% 
  group_by(clearance_status_d10, genus) %>% 
  mutate(median=(median(agg_rel_abund))) %>% #create a column of median values for each group
  ungroup() %>% 
  ggplot(aes(x=clearance_status_d10, y =agg_rel_abund, colour= clearance_status_d10))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median, ymin = median), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(aes(shape = subset), size=2, show.legend = TRUE, alpha = .4) +
  scale_colour_manual(name=NULL,
                      values=c("seagreen3", "mediumpurple"),
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
        strip.text = element_text(hjust = 0.5, size = 6.8, face = "italic"),
        axis.text.x = element_blank(),
        legend.position = "bottom")
save_plot(filename = paste0("results/figures/ml_d5_top10_genus.png"), top10_d5_model_taxa, base_height = 5, base_width = 8)  

#Make composite figure of ML results for 5dpi----
top_panel <- plot_grid(performance_otu, rf_feat_5dpi, labels = NULL, label_size = 12, nrow=1, rel_widths = c(1,2))
plot_grid(top_panel, top10_d5_model_taxa, labels = NULL, label_size = 12, ncol=1)+
  ggsave("results/figures/ml_summary_5dpi_genus.pdf", width=8, height=8)

