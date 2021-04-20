source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Define color scheme for waters, mocks & FMTs----
color_scheme <- c("#67a9cf", "#ef8a62", "7f5f1e") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("water", "mock", "FMT")
color_labels <- c( "Water", "Mock", "FMT")

#Read in diversity data for water & mock samples (+FMT)----
diversity_data <- read_tsv("data/water_test/peg3350.opti_mcc.groups.ave-std.summary") %>%
  filter(method == "ave") %>%
  select(group, sobs, shannon, invsimpson, coverage) %>%
  rename(unique_label = group) %>% #group is the same as unique_label in the metadata data frame
  mutate(sample_type = case_when(str_detect(unique_label, "water") ~ "water",
                                 str_detect(unique_label, "FMT") ~ "FMT",
                                 str_detect(unique_label, "mock") ~ "mock",
                                 TRUE ~ unique_label)) %>%
  mutate(day = NA, #add columns needed for plot_pcoa function
         group = sample_type) #add columns needed for plot_pcoa function

#Function to plot different alpha diversity metrics with the following arguments
#alpha_metric: how alpha metric of choice is listed in dataframe. Ex. sobs, shannon, etc.
#y_axis_label: how you want to label the alpha metric on the plot. Ex. "Shannon Diversity Index"
plot_alpha_metric <- function(alpha_metric, y_axis_label){
  diversity_data %>%
    mutate(group = fct_relevel(group, "water", "mock", "FMT")) %>% #Make sure order is correct
    group_by(group) %>%
    mutate(median = median({{ alpha_metric }})) %>% #Create column of median values for each group
    ungroup() %>%
    ggplot(aes(x=group, y = {{ alpha_metric }}, color = group))+
    geom_errorbar(aes(ymax= median, ymin= median, color = group), size = 1)+#Add line to show median of each point
    geom_jitter(size=2, show.legend = FALSE) +
    labs(title=NULL,
         x=NULL,
         y=y_axis_label)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels,
                        guide = "none")+
    scale_x_discrete(label = c("Water", "Mock", "FMT"))+
    theme_classic()+
    theme(legend.position = "bottom",
          text = element_text(size = 19),# Change font size for entire plot
          axis.title.y = element_text(size = 17))
}

#Shannon, inverse simpson and richness plots
shannon_plot <- plot_alpha_metric(shannon, "Shannon Diversity Index")
save_plot("exploratory/notebook/waters_mock_fmts_shannon.png", shannon_plot)
richness_plot <- plot_alpha_metric(sobs, "Number of Observed OTUs")
save_plot("exploratory/notebook/waters_mock_fmts_richness.png", richness_plot)

#Read in PCoA data for water & mock samples (+FMT)----
pcoa_data <- read_tsv("data/water_test/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(unique_label = group) %>% #group is the same as id in the metadata data frame
  mutate(sample_type = case_when(str_detect(unique_label, "water") ~ "water",
                                 str_detect(unique_label, "FMT") ~ "FMT",
                                 str_detect(unique_label, "mock") ~ "mock",
                                 TRUE ~ unique_label)) %>%
  mutate(day = NA, #add columns needed for plot_pcoa function
         group = sample_type) #add columns needed for plot_pcoa function
#Read in .loadings file to add percent variation represented by PCoA axis
axis_labels <- read_tsv("data/water_test/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- axis_labels %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- axis_labels %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

pcoa <- plot_pcoa(pcoa_data)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))
save_plot("exploratory/notebook/waters_mock_fmts_pcoa.png", pcoa, base_height = 5, base_width = 5)

#Create stand alone legend
group_legend <- pcoa_data %>%
  ggplot(aes(x = axis1, y = axis2, color = group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  geom_point()+ theme_classic()
group_legend <- get_legend(group_legend)
save_plot("exploratory/notebook/waters_mock_fmts_legend.png", group_legend, base_height = .8, base_width = .8)

#Read in shared and taxonomy files for water & mock samples (+FMT)----
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

#Examine Peptostreptococcaceae OTU 2018 in water, mocks----
#(Different OTU number because I analyzed these samples separately with mothur, otherwise the mock samples are filtered out)
pepto_otu <- agg_otu_data %>% 
  mutate(group = case_when(str_detect(unique_label, "water") ~ "water", #Create groups based on sample types
                                 str_detect(unique_label, "FMT") ~ "FMT",
                                 str_detect(unique_label, "mock") ~ "mock",
                                 TRUE ~ unique_label)) %>%
  filter(otu == "Peptostreptococcaceae (OTU 2018)") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% # 2,000 is 2 times the subsampling parameter of 1000
  mutate(group = fct_relevel(group, "water", "mock", "FMT")) %>% #Make sure order is correct
  group_by(group) %>%
  mutate(median = median(agg_rel_abund)) %>% #Create column of median values for each group
  ungroup() %>%
  ggplot(aes(x=group, y = agg_rel_abund, color = group))+
  geom_errorbar(aes(ymax= median, ymin= median, color = group), size = 1)+#Add line to show median of each point
  geom_jitter(size=2, show.legend = FALSE) +
  geom_hline(yintercept=1/1000, color="gray")+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels,
                      guide = "none")+
  scale_x_discrete(label = c("Water", "Mock", "FMT"))+
  theme_classic()+
  labs(title=NULL,
        x=NULL,
        y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme(legend.position = "bottom",
        text = element_text(size = 19),# Change font size for entire plot
        axis.title.y = element_text(size = 17))
ggsave("exploratory/notebook/pepto_otu_water_mock.png", pepto_otu)

#PCoA data for all samples together----
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


