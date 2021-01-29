source("code/utilities.R") #Loads libraries, reads in metadata, functions

#List of unique_labels corresponding to the stool samples for the tissue samples that were sequenced
stool_corresponding_to_tissues <- metadata %>% 
  filter(sample_type == "cecum") %>% #select one of the 3 types of tissue that were sequenced for each mouse
  separate(unique_label, into = c("tissue", "corresponding_stool_sample"), sep = 1) %>% 
  pull(corresponding_stool_sample)

#Read in PCoA for all samples----
all_pcoa_data <- read_tsv("data/process/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(unique_label = group) %>% #group is the same as id in the metadata data frame
  left_join(metadata, by= "unique_label")  #merge metadata and PCoA data frames


#Read in .loadings file to add percent variation represented by PCoA axis
all_axis_labels <- read_tsv("data/process/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
all_axis1 <- all_axis_labels %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
all_axis2 <- all_axis_labels %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

#Plot PCoA data with the 3 waters that had > 1000 sequences 
pcoa_water <- all_pcoa_data %>% 
  mutate(group = ifelse(grepl("water", unique_label), "water", "sample")) %>%  #Temporarily create just 2 groups. Waters and the rest of the samples
  ggplot(aes(x=axis1, y=axis2, color = group))+
  geom_point(size=2)+
  theme_classic()+
  theme(legend.position = c(.9, .2))+
  theme(legend.position = c(.9, .3))+
  xlim(-0.425, 0.65)+
  ylim(-0.525, 0.5)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))
save_plot(filename = paste0("exploratory/notebook/pcoa_samples_water.png"), pcoa_water, base_height = 5, base_width = 5)

#PCoA data of tissue samples collected from the experiment endpoint with their corresponding stool samples
pcoa_tissues_v_stool <- all_pcoa_data %>% 
  filter(sample_type == "water"| sample_type == "cecum"| sample_type == "distal_colon"| sample_type == "proximal_colon"|unique_label %in% stool_corresponding_to_tissues) %>% #select the tissue samples & the corresponding stool samples that were collected from the same mice at the same timepoint
  ggplot(aes(x=axis1, y=axis2, color = sample_type))+
  geom_point(size=2)+
  theme_classic()+
  theme(legend.position = c(.9, .3))+
  xlim(-0.425, 0.65)+
  ylim(-0.525, 0.5)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))
save_plot(filename = paste0("exploratory/notebook/pcoa_tissues_v_stools.png"), pcoa_tissues_v_stool, base_height = 5, base_width = 5)

#Read in alpha diversity metrics for all samples----
diversity_data <- read_tsv("data/process/peg3350.opti_mcc.groups.ave-std.summary") %>%
  filter(method == "ave") %>%
  select(group, sobs, shannon, invsimpson, coverage) %>%
  rename(unique_label = group) %>% #group is the same as unique_label in the metadata data frame
  left_join(metadata, by = "unique_label") #Match only the samples we have sequence data for

#Read in shared and taxonomy files for all samples----
# Import taxonomy into data frame and clean up taxonomy names
taxonomy <- read_tsv(file="data/process/peg3350.taxonomy") %>%
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
otu_data <- read_tsv("data/process/peg3350.opti_mcc.0.03.subsample.shared", col_types=cols(Group=col_character())) %>%
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
agg_genus_data <- agg_taxonomic_data(genus)

#Rename otus to match naming convention used previously and add a column that will work with ggtext package:
agg_otu <- agg_otu_data %>%
  mutate(key=otu) %>%
  group_by(key)
taxa_info <- read.delim('data/process/peg3350.taxonomy', header=T, sep='\t') %>%
  select(-Size) %>%
  mutate(key=OTU) %>%
  select(-OTU)
agg_otu_data <- inner_join(agg_otu, taxa_info, by="key") %>%
  ungroup() %>%
  mutate(key=str_to_upper(key)) %>%
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>%
  mutate(taxa=gsub("(.*)_.*","\\1",Taxonomy)) %>%
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>%
  mutate(taxa=str_replace_all(taxa, c("Clostridium_" = "Clostridium "))) %>% 
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

#Remove environmental variables that are no longer needed
rm(agg_otu, agg_taxa_data, duplicated_seq_samples, duplicates_to_drop, otu_data, peg3350.files,
     contaminated_notes, contaminated_samples, seq_files_missing_from_metadata, prep_notes)
