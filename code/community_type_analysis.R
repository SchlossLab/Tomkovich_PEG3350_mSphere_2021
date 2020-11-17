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
sample_best_community_fit <- read_tsv("data/process/peg3350.subsample.genus.genus.dmm.mix.design", col_names=c("unique_label", "best_fitting_partition"))
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
bacteria_in_clusters <- community_otus %>% 
  left_join(taxonomy, by = "otu") %>% #Join to taxonomy data frame by otu column in order to get the genus name for each OTU
  select(genus,ends_with("mean")) %>% 
  select(-P0.mean) %>% 
  slice_head(n=10) %>% 
  pivot_longer(-genus,names_to="cluster",values_to="relabund") %>% 
  mutate(cluster=str_replace(cluster,"P(\\d*).mean","\\1")) %>% 
  ggplot() + geom_tile(aes(x=genus,y=cluster,fill=relabund))+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic"))
  
#Examine the sample types that belong to each community type:

sample_type_clustering <- sample_best_community_fit %>% 
  left_join(metadata, by = "unique_label") %>% 
  group_by(sample_type) %>% 
  count(best_fitting_partition) %>% 
  ggplot()+
  geom_boxplot(aes(x=sample_type, y=n, 
                   color = best_fitting_partition))+
  theme_classic()

sample_type_clustering_v2 <- sample_best_community_fit %>% 
  left_join(metadata, by = "unique_label") %>% 
  group_by(sample_type) %>% 
  count(best_fitting_partition) %>% 
  ggplot()+
  geom_boxplot(aes(x=best_fitting_partition, y=n, 
                   color = sample_type))+
  coord_flip()+
  theme_classic()+
  geom_vline(xintercept = c((1:15) - 0.5 ), color = "grey")  # Add gray lines to clearly separate partitions
