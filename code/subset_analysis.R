source("code/utilities.R")

#Check for overlap in CFU/Weight figs----------
#inner join to check for overlap (no overlap present)
one_day_PEG_metadata %>% inner_join(five_day_PEG_metadata)
one_day_PEG_metadata %>% rbind(five_day_PEG_metadata) %>% inner_join(post_cdi_PEG_metadata)


#Split sequenced samples into 1 day PEG, 5 day PEG, and post-CDI subsets----------

#Full list of samples with sequence data that includes all timepoints
pcoa_data <- read_tsv("data/process/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>% 
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(unique_label = group) %>% #group is the same as unique_label in the metadata data frame
  right_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

#1 day PEG samples----
pcoa_1_day_PEG <- pcoa_data %>% 
  inner_join(one_day_PEG_metadata) %>% 
  pull(unique_label) %>%
  noquote() #Remove quotations from the characters

#Concatenate output and add - between each sample.
one_day_PEG_unique_labels <- (paste(pcoa_1_day_PEG, collapse = "-"))

#5 day PEG samples----
pcoa_5_day_PEG <- pcoa_data %>% 
  inner_join(five_day_PEG_metadata) %>% 
  pull(unique_label) %>%
  noquote() #Remove quotations from the characters

#Concatenate output and add - between each sample.
five_day_PEG_unique_labels <- (paste(pcoa_5_day_PEG, collapse = "-"))

#5 day PEG samples----
pcoa_post_CDI_PEG <- pcoa_data %>% 
  inner_join(post_cdi_PEG_metadata) %>% 
  pull(unique_label) %>%
  noquote() #Remove quotations from the characters

#Concatenate output and add - between each sample.
post_CDI_PEG_unique_labels <- (paste(pcoa_post_CDI_PEG, collapse = "-"))

