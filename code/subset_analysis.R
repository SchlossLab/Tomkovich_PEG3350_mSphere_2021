source("code/utilities.R")

#Check for overlap in CFU/Weight figs----------
#Create the same metadata subsets used in cfu_weight.R scripts for the 3 subsets----
one_day_PEG_metadata <- one_day_PEG_subset(metadata)
five_day_PEG_metadata <- five_day_PEG_subset(metadata) %>%
  filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon")) #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected
post_cdi_PEG_metadata <- post_cdi_PEG_subset(metadata) %>%
  filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon"))  #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected

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

#Create group argument to use in mothur dist.shared() function----
#Use the following function to remove spacing and quotation around each sample's unique_label
#df = dataframe representing a subset of samples to put into the group= option in mothur's dist.shared function
mothur_format <- function(df){
  df %>% 
    pull(unique_label) %>%
    noquote() %>%  #Remove quotations from the characters
    glue_collapse(sep = "-") #Separate all samples with a -
}

#1 day PEG samples----
pcoa_1_day_PEG <- pcoa_data %>% 
  inner_join(one_day_PEG_metadata)

#Concatenate output and add - between each sample.
one_day_PEG_unique_labels <- mothur_format(pcoa_1_day_PEG)

#5 day PEG samples with mock & 5 day PEG tissues added----
pcoa_5_day_PEG <- pcoa_data %>% 
  inner_join(five_day_PEG_metadata) %>% 
  add_row(tissue_5_days_PEG) %>% #Add tissues
  add_row(mock_stools_tissues) 

#Concatenate output and add - between each sample.
five_day_PEG_unique_labels <- mothur_format(pcoa_5_day_PEG)

#Split 5 day PEG into stool and tissue subsets:
five_day_PEG_stool_labels <- mothur_format(pcoa_5_day_PEG %>% filter(sample_type == "stool") %>% 
                                             filter(!group %in% c("WMN", "CN"))
                                           )
  
five_day_PEG_tissue_labels <-   mothur_format(pcoa_5_day_PEG %>% filter(!sample_type == "stool")%>% 
                                                filter(!group %in% c("WMN", "CN"))
                                              )

#Create subsets with mock samples
five_day_PEG_stool_mock_labels <- mothur_format(pcoa_5_day_PEG %>% filter(sample_type == "stool") %>% 
                                                  #For these PCoAs we just want to compare clindamycin +/- C. difficile and 5-day PEG +/- C. difficilie
                                                  filter(!group %in% c("WMC", "WMR")))

five_day_PEG_tissue_mock_labels <- mothur_format(pcoa_5_day_PEG %>% filter(!sample_type == "stool")%>% 
                                                   #For these PCoAs we just want to compare clindamycin +/- C. difficile and 5-day PEG +/- C. difficilie
                                                   filter(!group %in% c("WMC", "WMR"))) 

#Post-CDI PEG samples with FMT solution & Post-CDI PEG tissues added----
pcoa_post_CDI_PEG <- pcoa_data %>% 
  inner_join(post_cdi_PEG_metadata) %>% 
  add_row(FMT_subset) %>% 
  add_row(tissues_post_CDI_PEG)

#Concatenate output and add - between each sample.
post_CDI_PEG_unique_labels <- mothur_format(pcoa_post_CDI_PEG)

#Split post-CDI PEG into stool (+FMT) and tissue subsets:
post_CDI_PEG_stool_labels <- mothur_format(pcoa_data %>% 
                                             inner_join(post_cdi_PEG_metadata) %>% 
                                             add_row(FMT_subset))
post_CDI_PEG_tissue_labels <- mothur_format(tissues_post_CDI_PEG)
