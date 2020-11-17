source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Read in shared and taxonomy files for water & mock samples (+FMT & PBS gavage)----
# Import taxonomy into data frame and clean up taxonomy names
taxonomy <- read_tsv(file="data/water_test/peg3350.taxonomy") %>%
  rename_all(tolower) %>% #remove uppercase from column names
  # Split taxonomic information into separate columns for each taxonomic level
  mutate(taxonomy=str_replace_all(taxonomy, c("\\(\\d*\\)" = "", #drop digits with parentheses around them
                                              ';$' = "", #removes semi-colon at end of line
                                              'Bacteria_unclassified' = 'Unclassified',
                                              "Clostridium_" = "Clostridium ", #Remove underscores after Clostridium
                                              "_" = " ", #Removes all other underscores
                                              "unclassified" = "Unclassified"))) %>%
  # Separate taxonomic levels into separate columns according to semi-colon.
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';')
#Warning message: ~80 OTUs that are unknown, filled in with NAs

# Import otu_data for samples
otu_data <- read_tsv("data/water_test/peg3350.opti_mcc.0.03.subsample.shared", col_types=cols(Group=col_character())) %>%
  select(-label, -numOtus) %>%
  rename(unique_label = Group) %>% #group is the same as unique_label in the metadata data frame
  gather(-unique_label, key="otu", value="count") %>%
  mutate(rel_abund=count/1000) #Use 1000, because this is the subsampling parameter chosen.

#Merge otu_data to taxonomy data frame
agg_taxa_data <- inner_join(otu_data, taxonomy)

# Function to summarize relative abundance level for a given taxonomic level (ex. genus, family, phlyum, etc.)
agg_taxonomic_data <- function(taxonomic_level) {
  agg_taxa_data %>%
    group_by(unique_label, {{ taxonomic_level }}) %>% #Embracing treats the taxonomic_level argument as a column name
    summarize(agg_rel_abund=sum(rel_abund)) %>%
    # Merge relative abundance data to specifci taxonomic_level data
    left_join(., metadata, by = "unique_label") %>%
    ungroup()
}

# Relative abundance data at the otu level:
agg_otu_data <- agg_taxonomic_data(otu) 

#Rename otus to match naming convention used previously and add a column that will work with ggtext package:
agg_otu <- agg_otu_data %>%
  mutate(key=otu) %>%
  group_by(key)
taxa_info <- read.delim('data/water_test/peg3350.taxonomy', header=T, sep='\t') %>%
  select(-Size) %>%
  mutate(key=OTU) %>%
  select(-OTU)
agg_otu_data <- inner_join(agg_otu, taxa_info, by="key") %>%
  ungroup() %>%
  mutate(key=str_to_upper(key)) %>%
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>%
  mutate(taxa=gsub("(.*)_.*","\\1",Taxonomy)) %>%
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>%
  mutate(taxa=gsub(".*;","",taxa)) %>%
  mutate(taxa=gsub("(.*)_.*","\\1",taxa)) %>%
  mutate(taxa=gsub('[0-9]+', '', taxa)) %>%
  mutate(taxa=str_remove_all(taxa, "[(100)]")) %>%
  unite(key, taxa, key, sep=" (") %>%
  mutate(key = paste(key,")", sep="")) %>%
  select(-otu, -Taxonomy) %>%
  rename(otu=key) %>%
  mutate(otu=paste0(gsub('TU0*', 'TU ', otu))) %>%
  separate(otu, into = c("bactname", "OTUnumber"), sep = "\\ [(]", remove = FALSE) %>% #Add columns to separate bacteria name from OTU number to utilize ggtext so that only bacteria name is italicized
  mutate(otu_name = glue("*{bactname}* ({OTUnumber}")) #Markdown notation so that only bacteria name is italicized


# Examination of the OTUs in the water controls with > 1000 sequences----
#Figure out the top 10 OTUs across the water controls:
top_10_water_otus <- agg_otu_data %>% 
  filter(str_detect(unique_label, "water")) %>% 
  mutate(sample_type = "water") %>% 
  group_by(sample_type, otu) %>% 
  summarise(mean_rel_abund = mean(agg_rel_abund)) %>% 
  top_n(10, mean_rel_abund) %>% 
  pull(otu)
#Plot the top 10 OTUs across the 3 water controls:
water_otus <- agg_otu_data %>% 
  filter(str_detect(unique_label, "water")) %>% #select only the water controls
  filter(otu %in% top_10_water_otus) %>% #Select the top 10 water OTUs
  ggplot(aes(fill = otu, y = agg_rel_abund, x = unique_label))+
  geom_bar(position = "fill", stat="identity")+
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance")+
  coord_flip()+
  theme_classic()

#Figure out the top 10 OTUs across the mock controls:
top_10_mock_otus <- agg_otu_data %>% 
  filter(str_detect(unique_label, "mock")) %>% 
  mutate(sample_type = "mock") %>% 
  group_by(sample_type, otu) %>% 
  summarise(mean_rel_abund = mean(agg_rel_abund)) %>% 
  top_n(10, mean_rel_abund) %>% 
  pull(otu)
#Plot the top 10 OTUs across the mock controls:
mock_otus <- agg_otu_data %>% 
  filter(str_detect(unique_label, "mock")) %>% #select only the mock controls
  filter(otu %in% top_10_mock_otus) %>% #Select the top 10 mock OTUs
  ggplot(aes(fill = otu, y = agg_rel_abund, x = unique_label))+
  geom_bar(position = "fill", stat="identity")+
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance")+
  coord_flip()+
  theme_classic()

#All samples with top 10 mock otus:
all_samples_mock_otus <- agg_otu_data %>% 
  filter(otu %in% top_10_mock_otus) %>% #Select the top 10 water OTUs
  ggplot(aes(fill = otu, y = agg_rel_abund, x = unique_label))+
  geom_bar(position = "fill", stat="identity")+
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance")+
  coord_flip()+
  theme_classic()

# of sequences in the mock and water controls for plate 14, 15, 16, and 17
#Plate 14: Mock 24892 Water 4387
#Plate 15: Mock 28474 Water 253
#Plate 16: Mock 174716 Water 5122
#Plate 17: Mock 41207 Water 1683

#All samples top 30 otus
top_30_otus <- agg_otu_data %>% 
  group_by(otu) %>% 
  summarise(mean_rel_abund = mean(agg_rel_abund)) %>% 
  top_n(30, mean_rel_abund) %>% 
  pull(otu)

#All samples with top 30 mock otus:
all_samples_otus <- agg_otu_data %>% 
  filter(otu %in% top_30_otus) %>% 
  ggplot(aes(fill = otu, y = agg_rel_abund, x = unique_label))+
  geom_bar(position = "fill", stat="identity")+
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance")+
  coord_flip()+
  theme_classic()
save_plot(filename = "exploratory/notebook/waters_mock_fmts_top_30_otus.png", all_samples_otus, base_height = 4.5, base_width = 8.5, base_aspect_ratio = 2)


