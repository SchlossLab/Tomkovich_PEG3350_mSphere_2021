source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Prepare input files for community type analysis with mothur----
#Create .shared and .taxonomy files at the genus level to use in Dirichlet Multinomial Mixture analysis
#Use these as input for mothur's get.communitytype function. See code/community_type.batch
#Shared file:
shared <- read.delim('data/process/peg3350.opti_mcc.0.03.subsample.shared', header=T, sep='\t') %>%
  select(-label, -numOtus) %>% 
  gather(-Group, key=OTU, value=count)

#Read in taxonomy and select genus level:
taxonomy <- read_tsv(file="data/process/peg3350.taxonomy") %>% 
  select(-Size) %>%
  mutate(Taxonomy=str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>%
  mutate(Taxonomy=str_replace_all(Taxonomy, ";$", "")) %>%
  separate(Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';') %>%
  select(OTU, "genus") %>%
  rename(taxon = genus)

unique_taxonomy <- taxonomy %>%
  select(taxon) %>%
  unique() %>%
  mutate(otu = paste0("Otu", str_pad(1:nrow(.), width=nchar(nrow(.)), pad="0")))

#Join genus level taxonomy to shared to create shared file at the genus level:
genus_shared <- inner_join(shared, taxonomy, by="OTU") %>%
  group_by(taxon, Group) %>%
  summarize(count = sum(count)) %>%
  ungroup() %>%
  inner_join(., unique_taxonomy) %>%
  select(-taxon) %>%
  spread(otu, count) %>%
  mutate(label="genus", numOtus=ncol(.)-1) %>%
  select(label, Group, numOtus, everything())
write_tsv(genus_shared, path = "data/process/peg3350.subsample.genus.shared")

select(genus_shared, -label, -numOtus) %>%
  gather(otu, count, -Group) %>%
  group_by(otu) %>%
  summarize(count=sum(count)) %>%
  inner_join(., unique_taxonomy) %>%
  rename("OTU"="otu", "Size"="count", "Taxonomy"="taxon") %>%
  write_tsv(path ="data/process/peg3350.genus.taxonomy")

#Visualize get.communitytype analysis results----
#Read in dmm fit data to evaluate community type fit with the laplace value depending on the number of community types
dmm_fit <- read_tsv("data/process/peg3350.subsample.genus.genus.dmm.mix.fit")

laplace_plot <- dmm_fit %>% 
  ggplot()+
  geom_line(aes(x = K, y = Laplace))+
  theme_classic()
#Save results
save_plot(filename = "exploratory/notebook/motility_community_type_laplace.png", laplace_plot)

#Read in the rest of the output files from mothur's get.communitytype function:
#Read in data that indicate the bacteria (at the genus level) membership of the 15 different  community types
community_otus <- read_tsv("data/process/peg3350.subsample.genus.genus.dmm.mix.summary") %>% 
  rename("otu" = "OTU") #rename to match otu column in taxonomy
community_parameters <- read_tsv("data/process/peg3350.subsample.genus.genus.dmm.mix.parameters")
sample_community_membership <- read_tsv("data/process/peg3350.subsample.genus.genus.dmm.15.mix.posterior")
sample_best_community_fit <- read_tsv("data/process/peg3350.subsample.genus.genus.dmm.mix.design", col_names=c("unique_label", "best_fitting_community")) %>% 
  left_join(metadata, by = "unique_label") %>%  #Join to metadata
  mutate(community=str_replace(best_fitting_community,"Partition_(\\d*)","\\1")) %>% #Add a column with just the community number
  mutate(community= as.numeric(community)) #Transform community variable type from character to numeric
sample_community_rel_abund <- read_tsv("data/process/peg3350.subsample.genus.genus.dmm.15.mix.relabund")

## Import genus-level taxonomy into data frame and clean up names
taxonomy <- read_tsv(file="data/process/peg3350.genus.taxonomy") %>%
  rename_all(tolower) %>% #remove uppercase from column names
  rename(genus = taxonomy) %>% #Rename taxonomy to genus
  # Clean up genus names  
  mutate(genus=str_replace_all(genus, c('Bacteria_unclassified' = 'Unclassified',
                                        "Clostridium_" = "Clostridium ", #Remove underscores after Clostridium
                                        "_" = " ", #Removes all other underscores
                                        "unclassified" = "unclassified")))

#Examine the bacteria that make up each community type:
bacteria_in_community <- community_otus %>% 
  left_join(taxonomy, by = "otu") %>% #Join to taxonomy data frame by otu column in order to get the genus name for each OTU
  select(genus,ends_with("mean")) %>% 
  select(-P0.mean) %>% 
  slice_head(n=16) %>% #Top 16 includes Peptostreptococcaceae, which is likely C. difficile
  pivot_longer(-genus,names_to="community",values_to="relabund") %>% 
  mutate(community=str_replace(community,"P(\\d*).mean","\\1")) %>% 
  mutate(community= as.numeric(community)) %>% #Transform community variable type from character to numeric
  ggplot() + geom_tile(aes(x=genus,y=community,fill=relabund))+
  coord_flip()+
  scale_y_continuous(breaks = c(1:15))+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic"))+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Relative \nAbundance")+
  theme(axis.title.y = element_blank()) #Get rid of y axis title
  
#Examine the sample types and groups that belong to each community type:
sample_type_clustering <- sample_best_community_fit %>% 
  group_by(sample_type) %>% 
  count(community) %>% 
  ggplot()+
  geom_jitter(aes(x=community, y=n, 
                   color = sample_type))+
  scale_x_continuous(breaks = c(1:15))+
  coord_flip()+
  theme_classic()+
  geom_vline(xintercept = c((1:15) - 0.5 ), color = "grey")  # Add gray lines to clearly separate partitions

#Obtain total number of samples per cluster
sample_best_community_fit %>% 
  group_by(community) %>% 
  summarize(cluster_total = n())

percent_sample_type <- sample_best_community_fit %>% 
  group_by(community, sample_type) %>% 
  summarize(sample_type_community_total = n()) %>%
  #Make a new variable % group per community type based on total community numbers
  mutate(percent_community = case_when(community == "1" ~ (sample_type_community_total/74)*100,
                                     community == "2" ~ (sample_type_community_total/80)*100,
                                     community == "3" ~ (sample_type_community_total/63)*100,
                                     community == "4" ~ (sample_type_community_total/112)*100,
                                     community == "5" ~ (sample_type_community_total/48)*100,
                                     community == "6" ~ (sample_type_community_total/37)*100,
                                     community == "7" ~ (sample_type_community_total/189)*100,
                                     community == "8" ~ (sample_type_community_total/120)*100,
                                     community == "9" ~ (sample_type_community_total/167)*100,
                                     community == "10" ~ (sample_type_community_total/42)*100,
                                     community == "11" ~ (sample_type_community_total/82)*100,
                                     community == "12" ~ (sample_type_community_total/100)*100,
                                     community == "13" ~ (sample_type_community_total/33)*100,
                                     community == "14" ~ (sample_type_community_total/143)*100,
                                     community == "15" ~ (sample_type_community_total/94)*100,
                                     TRUE ~ 0)) %>% #No samples should fall into this category
  mutate(sample_type = factor(sample_type, levels = unique(as.factor(sample_type)))) %>% #Transform into factor variable
  mutate(sample_type = fct_relevel(sample_type, "water", "FMT", "stool", "distal_colon", "proximal_colon", "cecum")) %>% #Specify the order of the groups
  ggplot()+
  geom_tile(aes(x=community, y=sample_type, fill=percent_community))+
  scale_x_continuous(breaks = c(1:15))+
  scale_y_discrete(label = c("Water", "FMT", "Stool", "Distal colon", "Proximal colon", "Cecum"))+
  theme_classic()+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "% Community")+
  theme(axis.title.y = element_blank(), #Get rid of y axis title
        axis.title.x = element_blank(), #Get rid of x axis title, text, and ticks. Will combine with bacteria in community
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

percent_group <- sample_best_community_fit %>% 
  group_by(community, group) %>% 
  summarize(group_community_total = n()) %>%
  #Make a new variable % group per community type based on total community numbers
  mutate(percent_community = case_when(community == "1" ~ (group_community_total/74)*100,
                                       community == "2" ~ (group_community_total/80)*100,
                                       community == "3" ~ (group_community_total/63)*100,
                                       community == "4" ~ (group_community_total/112)*100,
                                       community == "5" ~ (group_community_total/48)*100,
                                       community == "6" ~ (group_community_total/37)*100,
                                       community == "7" ~ (group_community_total/189)*100,
                                       community == "8" ~ (group_community_total/120)*100,
                                       community == "9" ~ (group_community_total/167)*100,
                                       community == "10" ~ (group_community_total/42)*100,
                                       community == "11" ~ (group_community_total/82)*100,
                                       community == "12" ~ (group_community_total/100)*100,
                                       community == "13" ~ (group_community_total/33)*100,
                                       community == "14" ~ (group_community_total/143)*100,
                                       community == "15" ~ (group_community_total/94)*100,
                                       TRUE ~ 0)) %>% #No samples should fall into this category
  mutate(group = fct_relevel(group, "water", "FMT", "CN", "C", "FRM", "RM", "CWM", "1RM1", "M1", "WMR", "WMC", "WMN", "WM")) %>% #Specify the order of the groups
  ggplot()+
  geom_tile(aes(x=community, y=group, fill=percent_community))+
  scale_x_continuous(breaks = c(1:15))+
  scale_y_discrete(label = c("Water", "FMT", "Clind. without infection", "Clind.", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350",
                             "Clind. + 1-day PEG 3350", "1-day PEG 3350 + 1-day recovery", "1-day PEG 3350", "5-day PEG 3350 + 10-day recovery", "5-day PEG 3350 + Clind.", 
                             "5-day PEG 3350 without infection", "5-day PEG 3350"))+ #Descriptive group names that match the rest of the plots
  theme_classic()+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "% Community")+
  theme(axis.title.y = element_blank(), #Get rid of y axis title
        axis.title.x = element_blank(), #Get rid of x axis title, text, and ticks. Will combine with bacteria in community
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())  

percent_day <- sample_best_community_fit %>% 
  group_by(community, day) %>% 
  summarize(day_community_total = n()) %>%
  #Make a new variable % day per community type based on total community numbers
  mutate(percent_community = case_when(community == "1" ~ (day_community_total/74)*100,
                                       community == "2" ~ (day_community_total/80)*100,
                                       community == "3" ~ (day_community_total/63)*100,
                                       community == "4" ~ (day_community_total/112)*100,
                                       community == "5" ~ (day_community_total/48)*100,
                                       community == "6" ~ (day_community_total/37)*100,
                                       community == "7" ~ (day_community_total/189)*100,
                                       community == "8" ~ (day_community_total/120)*100,
                                       community == "9" ~ (day_community_total/167)*100,
                                       community == "10" ~ (day_community_total/42)*100,
                                       community == "11" ~ (day_community_total/82)*100,
                                       community == "12" ~ (day_community_total/100)*100,
                                       community == "13" ~ (day_community_total/33)*100,
                                       community == "14" ~ (day_community_total/143)*100,
                                       community == "15" ~ (day_community_total/94)*100,
                                       TRUE ~ 0)) %>% #No samples should fall into this category
  mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day into factor variable
  mutate(day = fct_relevel(day, NA, "30", "25", "20", "15", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1", "0",
                             "-1", "-2", "-4", "-5", "-10", "-11", "-15")) %>% #Specify the order of the groups
  ggplot()+
  geom_tile(aes(x=community, y=day, fill=percent_community))+
  scale_x_continuous(breaks = c(1:15))+
  theme_classic()+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "% Community")+
  theme(axis.title.y = element_blank(), #Get rid of y axis title
        axis.title.x = element_blank(), #Get rid of x axis title, text, and ticks. Will combine with bacteria in community
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#3 waters all belong to community type 13
# 3 FMTs belong to community types 1, 8, and 11

#Combine plots of percent of each group or sample type that belong to each cluster with the
#plot of the bacteria relative abundances within each cluster
#Align heat maps
combined_heatmaps_sample_type <- align_plots(percent_sample_type, bacteria_in_community, align = 'v', axis = "b")
combined_heatmaps_group <- align_plots(percent_group, bacteria_in_community, align = 'v', axis = "b")
combined_heatmaps_day <- align_plots(percent_day, bacteria_in_community, align = 'v', axis = "b")
#Plot and save the aligned heat maps
plot_grid(combined_heatmaps_sample_type[[1]], combined_heatmaps_sample_type[[2]], rel_heights = c(1, 3), ncol = 1)+
  ggsave("results/figures/community_types_sample_type.png", height = 6.0, width = 8.5)
plot_grid(combined_heatmaps_group[[1]], combined_heatmaps_group[[2]], rel_heights = c(1, 1.2), ncol = 1)+
  ggsave("results/figures/community_types_group.png", height = 8.0, width = 8.5)
plot_grid(combined_heatmaps_day[[1]], combined_heatmaps_day[[2]], rel_heights = c(1.5, 1), ncol = 1)+
  ggsave("results/figures/community_types_day.png", height = 10.5, width = 8.5)

#Community types across groups, facetted by experimental milestones: Baseline, Post-treatment, Post-treatment-----
#First transform timepoints to milestones
group_time_fit <- sample_best_community_fit %>% 
  filter(!group %in% c("FMT", "water")) %>% #Remove 3 water & 3 FMT samples for this analysis
  #Transform day variable into an integer variable to help with creating timepoint milestones
  mutate(timepoint_status = case_when(group == "WMR" & day == -15 ~ "baseline",
                                      group == "WMR" & day %in% c(-10, -5, -4, -2, -1, 0) ~ "post-treatment",
                                      group %in% c("C", "CWM", "RM", "FRM") & day < 0 ~"baseline",
                                      group %in% c("C", "CWM", "RM", "FRM") & day == 0 ~"post-treatment",
                                      group %in% c("WM", "WMC") & day == -5 ~"baseline",
                                      group %in% c("WM", "WMC") & day == -1 ~ "treatment",
                                      group %in% c("WM", "WMC") & day == 0 ~ "post-treatment",
                                      group == "1RM1" & day == -2 ~ "baseline",
                                      group == "1RM1" & day == 0 ~ "post-treatment",
                                      group == "M1" & day %in% c(-1, -11) ~ "baseline",
                                      group == "M1" & day ==0 ~ "post-treatment",
                                      group %in% c("WMN", "CN") & day == -5 ~"baseline",
                                      group == "CN" & day == -1 ~"baseline",
                                      group == "WMN" & day == -1 ~ "treatment",
                                      group == "WMN" & day %in% c(0, 4, 6, 30) ~ "post-treatment",
                                      group == "CN" & day %in% c(0, 4) ~"post-treatment",
                                 TRUE ~ "post-infection")) %>% 
  select(exp_num, group, day, timepoint_status, community)

#Make heatmap of group communitty fit, faceted by timepoint status
#Make labels for facets
facet_labels <- c("Baseline", "Treatment", "Post-treatment", "Post-infection")
names(facet_labels) <- c("baseline", "treatment", "post-treatment", "post-infection")
group_time_plot <- group_time_fit %>% 
  mutate(timepoint_status = fct_relevel(timepoint_status, "baseline", "treatment", "post-treatment", "post-infection")) %>%   
  group_by(community, group, timepoint_status) %>% 
  summarize(group_community_total = n()) %>%
  #Make a new variable % group per community type based on total community numbers
  mutate(percent_community = case_when(community == "1" ~ (group_community_total/74)*100,
                                       community == "2" ~ (group_community_total/80)*100,
                                       community == "3" ~ (group_community_total/63)*100,
                                       community == "4" ~ (group_community_total/112)*100,
                                       community == "5" ~ (group_community_total/48)*100,
                                       community == "6" ~ (group_community_total/37)*100,
                                       community == "7" ~ (group_community_total/189)*100,
                                       community == "8" ~ (group_community_total/120)*100,
                                       community == "9" ~ (group_community_total/167)*100,
                                       community == "10" ~ (group_community_total/42)*100,
                                       community == "11" ~ (group_community_total/82)*100,
                                       community == "12" ~ (group_community_total/100)*100,
                                       community == "13" ~ (group_community_total/33)*100,
                                       community == "14" ~ (group_community_total/143)*100,
                                       community == "15" ~ (group_community_total/94)*100,
                                       TRUE ~ 0)) %>% #No samples should fall into this category
  mutate(group = fct_relevel(group, "CN", "C", "FRM", "RM", "CWM", "1RM1", "M1", "WMR", "WMC", "WMN", "WM")) %>% #Specify the order of the groups
  ggplot()+
  geom_tile(aes(x=community, y=group, fill=percent_community))+
  facet_wrap(~timepoint_status, labeller = labeller(timepoint_status = facet_labels), nrow = 1)+
  scale_x_continuous(breaks = c(1:15))+
  scale_y_discrete(label = c("Clind. without infection", "Clind.", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350",
                             "Clind. + 1-day PEG 3350", "1-day PEG 3350 + 1-day recovery", "1-day PEG 3350", "5-day PEG 3350 + 10-day recovery", "5-day PEG 3350 + Clind.", 
                             "5-day PEG 3350 without infection", "5-day PEG 3350"))+ #Descriptive group names that match the rest of the plots
  theme_classic()+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "% Community")+
  theme(axis.title.y = element_blank(), #Get rid of y axis title
        plot.background=element_blank(),
        #remove plot border
        panel.border=element_blank(),
        #set thickness of axis ticks
        axis.ticks=element_line(size=0.4),
        strip.background = element_blank(), #get rid of box around facet_wrap labels
        axis.title.x = element_blank()) #Get rid of x axis title, text, and ticks. Will combine with bacteria in community
save_plot(filename = "results/figures/community_types_group_timepoint_status.png", group_time_plot, base_height = 5, base_width = 16)
      
#Plot just post-infection timepoint with bacteria in each community
post_infection_plot <- group_time_fit %>% 
  mutate(timepoint_status = fct_relevel(timepoint_status, "baseline", "treatment", "post-treatment", "post-infection")) %>%   
  filter(timepoint_status == "post-infection") %>% 
  group_by(community, group, timepoint_status) %>% 
  summarize(group_community_total = n()) %>%
  #Make a new variable % group per community type based on total community numbers
  mutate(percent_community = case_when(community == "1" ~ (group_community_total/74)*100,
                                       community == "2" ~ (group_community_total/80)*100,
                                       community == "3" ~ (group_community_total/63)*100,
                                       community == "4" ~ (group_community_total/112)*100,
                                       community == "5" ~ (group_community_total/48)*100,
                                       community == "6" ~ (group_community_total/37)*100,
                                       community == "7" ~ (group_community_total/189)*100,
                                       community == "8" ~ (group_community_total/120)*100,
                                       community == "9" ~ (group_community_total/167)*100,
                                       community == "10" ~ (group_community_total/42)*100,
                                       community == "11" ~ (group_community_total/82)*100,
                                       community == "12" ~ (group_community_total/100)*100,
                                       community == "13" ~ (group_community_total/33)*100,
                                       community == "14" ~ (group_community_total/143)*100,
                                       community == "15" ~ (group_community_total/94)*100,
                                       TRUE ~ 0)) %>% #No samples should fall into this category
  mutate(group = fct_relevel(group, "C", "FRM", "RM", "CWM", "1RM1", "M1", "WMR", "WMC", "WM")) %>% #Specify the order of the groups
  ggplot()+
  geom_tile(aes(x=community, y=group, fill=percent_community))+
  facet_wrap(~timepoint_status, labeller = labeller(timepoint_status = facet_labels), nrow = 1)+
  scale_x_continuous(breaks = c(1:15))+
  scale_y_discrete(label = c("Clind.", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350",
                             "Clind. + 1-day PEG 3350", "1-day PEG 3350 + 1-day recovery", "1-day PEG 3350", "5-day PEG 3350 + 10-day recovery", "5-day PEG 3350 + Clind.", 
                             "5-day PEG 3350"))+ #Descriptive group names that match the rest of the plots
  theme_classic()+
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "% Community")+
  theme(axis.title.y = element_blank(), #Get rid of y axis title
        plot.background=element_blank(),
        strip.background = element_blank(), #get rid of box around facet_wrap labels
        axis.title.x = element_blank(), #Get rid of x axis title, text, and ticks. Will combine with bacteria in community
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
#Combine post-infection with bacteria in community types
#Combine plots of percent of each group or sample type that belong to each cluster with the
#plot of the bacteria relative abundances within each cluster
#Align heat maps
combined_heatmaps_group_pi <- align_plots(post_infection_plot, bacteria_in_community, align = 'v', axis = "b")
#Plot and save the aligned heat maps
plot_grid(combined_heatmaps_group_pi[[1]], combined_heatmaps_group_pi[[2]], rel_heights = c(1, 1.2), ncol = 1)+
  ggsave("results/figures/community_types_group_postinfection.png", height = 8.0, width = 8.5)


