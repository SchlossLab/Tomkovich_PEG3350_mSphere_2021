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
#Read in OTU level model performance:
otu_results <- read_perf_results("results/performance_results.csv", 
                              "otu") 

#Examine AUC results for each machine learning method----
auc_results <- otu_results %>% 
  group_by(method) %>% 
  summarize(median = median(AUC)) %>% 
  arrange(desc(median))

#Plot performance for all methods with OTU input data----
performance_otu <- otu_results %>% 
  filter(taxa_level == "otu") %>% 
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
#        axis.text.x = element_blank(), #Remove x axis text (all models are at the OTU level),
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
  taxa_info <- read.delim('data/process/peg3350.taxonomy', header=T, sep='\t') %>%
    select(-Size) %>%
    mutate(names=OTU) %>%
    select(-OTU)
  final_feat_imp <- inner_join(feat_imp, taxa_info, by="names") %>%
    ungroup() %>%
    mutate(names=str_to_upper(names)) %>%
    mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>%
    mutate(taxa=gsub("(.*)_.*","\\1",Taxonomy)) %>%
    mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>%
    mutate(taxa=str_replace_all(taxa, c("Clostridium_" = "Clostridium "))) %>% 
    mutate(taxa=gsub(".*;","",taxa)) %>%
    mutate(taxa=gsub("(.*)_.*","\\1",taxa)) %>%
    mutate(taxa=gsub('[0-9]+', '', taxa)) %>%
    mutate(taxa=str_remove_all(taxa, "[(100)]")) %>%
    unite(names, taxa, names, sep=" (") %>%
    mutate(names = paste(names,")", sep="")) %>%
    select(-Taxonomy) %>%
    rename(otu=names) %>%
    mutate(otu=paste0(gsub('TU0*', 'TU ', otu))) %>%
    separate(otu, into = c("bactname", "OTUnumber"), sep = "\\ [(]", remove = FALSE) %>% #Add columns to separate bacteria name from OTU number to utilize ggtext so that only bacteria name is italicized
    mutate(otu_name = glue("*{bactname}* ({OTUnumber}")) #Markdown notation so that only bacteria name is italicized
  return(final_feat_imp)
}

rf_feat <- read_feat_imp("results/combined_feature-importance_rf.csv") #%>% 
  #Transform bactname variable into factor 
  #mutate(bactname = factor(bactname, levels = unique(as.factor(bactname))))

#Function to get the top 20 features that have the largest impact on AUROC
#df = dataframe of feature importances for the 100 seeds
top_20 <- function(df){
  data_first_20 <- df %>% 
    group_by(otu, bactname) %>% 
    summarize(median = median(perf_metric_diff)) %>% #Get the median performance metric diff. for each feature
    arrange(desc(median)) %>% #Arrange from largest median to smallest
    head(20)
  
  return(data_first_20)
}

#Top 20 OTUs for each input dataset with the random forest model
rf_top_feat <- top_20(rf_feat) %>% pull(otu)

#Create color schemes based on OTUs with the same taxa name----
color_scheme_df <-  top_20(rf_feat) %>% 
  group_by(bactname) %>% 
  mutate(duplicate = n()) %>% #count number of features with same bacteria name
  #Create color scheme based on number of times taxa shows up in top 20
  mutate(color = case_when(duplicate == 1 ~ "#bababa",
                           duplicate == 2 ~ "#404040",
                           duplicate == 3 ~ "blue",
                           duplicate == 11 ~ "#ca0020",
                           TRUE ~ "yellow")) #Should catch all cases with above, update if any are yellow
color_scheme <- color_scheme_df$color
color_bact <- color_scheme_df$bactname 
legend_bact <- color_bact

#Function to filter to top OTUs for each pairwise comparison & plot results----
#df = dataframes of feature importances for all seeds
#top_otus = dataframes of top otus
#comp_name = name of comparison to title the plot (in quotes)
plot_feat_imp <- function(df, top_otus){
  df %>% 
    filter(otu %in% top_otus) %>% 
    ggplot(aes(fct_reorder(otu_name, -perf_metric_diff, .desc = TRUE), perf_metric_diff, color = bactname))+
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
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          strip.background = element_blank(),
          legend.position = "none")   
}

#Plot feature importances for the top OTUs for each comparison----
plot_feat_imp(rf_feat, rf_top_feat)+
  ggsave("results/figures/ml_top_features_otu.png", height = 5, width = 8)
