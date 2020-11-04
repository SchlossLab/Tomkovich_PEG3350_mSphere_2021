source("code/utilities.R")

#Create the same metadata subsets used in cfu_weight.R scripts for the 3 subsets----
one_day_PEG_metadata <- one_day_PEG_subset(metadata)
five_day_PEG_metadata <- five_day_PEG_subset(metadata) %>%
  filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon")) #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected
post_cdi_PEG_metadata <- post_cdi_PEG_subset(metadata) %>%
  filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon"))  #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected

#Check what samples are not accounted for in our 3 subsets of C. diff cfu and weight data
unaccounted <- metadata %>% 
  filter(!unique_label %in% c(one_day_PEG_metadata$unique_label, five_day_PEG_metadata$unique_label, post_cdi_PEG_metadata$unique_label))
#Confirm unaccounted are primarily tissue samples and stool samples from mock challenged mice
#Unaccounted samples from tissue samples
unaccounted %>% filter(sample_type %in% c("cecum", "distal_colon", "proximal_colon")) %>% 
  group_by(sample_type) %>% 
  count() #252 samples are tissues from mice that already have a corresponding stool sample 
#We exclude tissues from our cfu and weight analysis because all 84 of these mice are already counted with their corresponding stool sample
#Unaccounted samples that are stool samples from mock challenged mice
unaccounted %>% filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon")) %>% 
  group_by(group) %>% 
  count() #292 samples that are in CN or WMN group (mock challenged clindamycin or 5-day PEG treated mice)

