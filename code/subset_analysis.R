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
  left_join(metadata, by= "unique_label") %>% #merge metadata and PCoA data frames 
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

#16S samples not accounted for in 3 subsets----
other_subset <- pcoa_data %>% 
  filter(!unique_label %in% c(one_day_PEG_metadata$unique_label, five_day_PEG_metadata$unique_label, post_cdi_PEG_metadata$unique_label))
#308 samples left out of current subsets. Most of these are tissues which were automatically filtered out of metadata.
# 3 FMT samples that should be grouped into post-CDI PEG subset.
# 3 water controls (likely contaminated do to low biomass tissue sequencing)
#Revise how we create the main 3 metadata subsets. Currently, all tissue samples are automatically excluded.

#Tissue samples to include in 5-days PEG subset: 
tissue_5_days_PEG <- other_subset %>% 
  filter(group == "C" & exp_num %in% c("M3","M4", "M5", "M8")| #Only use C mice from these experiments. Allocated groups to figures based on paper outline.
           group == "WM" & exp_num %in% c("M3","M4", "M5", "M8")|
           group == "WMC" & exp_num %in% c("M3","M4")|
           group == "WMR" & exp_num %in% c("M5","M6"))
  
#Stool & tissue samples for mock challenged mice to include in 5-days PEG subset:
mock_stools_tissues <- other_subset %>% 
  filter(group == "CN" | group == "WMN")

#Tissue samples to include in post_CDI_PEG subset:
tissues_post_CDI_PEG <- other_subset %>% 
  filter(group == "C" & exp_num %in% c("M7","M9")| #Only use C mice from these experiments. Allocated groups to figures based on paper outline.
           group == "CWM" & exp_num %in% c("M6","M7", "M9")|
           group == "RM" & exp_num %in% c("M7","M9")|
           group == "FRM" & exp_num %in% c("M9")) 

#FMT samples to include in post_CDI_PEG subset:
FMT_subset <-  other_subset %>% 
  filter(str_detect(unique_label, "FMT"))

#Unaccounted samples after creating 
unaccounted <- other_subset %>% 
  filter(!unique_label %in% c(tissue_5_days_PEG$unique_label, mock_stools_tissues$unique_label, 
                              tissues_post_CDI_PEG$unique_label, FMT_subset$unique_label))
#Only the 3 water controls remain after creating the 4 additional subsets 

#1 day PEG samples----
pcoa_1_day_PEG <- pcoa_data %>% 
  inner_join(one_day_PEG_metadata) %>% 
  pull(unique_label) %>%
  noquote() #Remove quotations from the characters

#Concatenate output and add - between each sample.
one_day_PEG_unique_labels <- (paste(pcoa_1_day_PEG, collapse = "-"))

#5 day PEG samples with mock & 5 day PEG tissues added----
pcoa_5_day_PEG <- pcoa_data %>% 
  inner_join(five_day_PEG_metadata) %>% 
  add_row(tissue_5_days_PEG) %>% #Add tissues
  add_row(mock_stools_tissues) %>% #Add samples from mock challenged mice
  pull(unique_label) %>%
  noquote() #Remove quotations from the characters

#Concatenate output and add - between each sample.
five_day_PEG_unique_labels <- (paste(pcoa_5_day_PEG, collapse = "-"))

#Post-CDI PEG samples with FMT solution & Post-CDI PEG tissues added----
pcoa_post_CDI_PEG <- pcoa_data %>% 
  inner_join(post_cdi_PEG_metadata) %>% 
  add_row(FMT_subset) %>% 
  add_row(tissues_post_CDI_PEG) %>% 
  pull(unique_label) %>%
  noquote() #Remove quotations from the characters

#Concatenate output and add - between each sample.
post_CDI_PEG_unique_labels <- (paste(pcoa_post_CDI_PEG, collapse = "-"))

