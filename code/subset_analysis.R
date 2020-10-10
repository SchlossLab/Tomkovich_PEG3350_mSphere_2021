source("code/utilities.R")

#Check for overlap in CFU/Weight figs----------
#1 day PEG Subset
one_day_PEG_metadata <- metadata %>%
  filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon")) %>% #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected
  filter(group == "C" & exp_num %in% c("M6")| #Only use C mice from this experiments. Allocated groups to figures based on paper outline.
           group == "1RM1" & exp_num %in% c("M6R")| #Had to differentiate experiment 6 from 6R in the metadata to create unique_mouse_id that wouldn't overlap for the M1 & 1RM1 mice that are both labeled with mouse_ids that are #s1-6
           group == "M1" & exp_num %in% c("M6"))%>%
  mutate(group=factor(group, levels=c("C", "1RM1", "M1")))



#5 days PEG Subset
five_day_PEG_metadata <- metadata %>%
  filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon")) %>% #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected
  filter(group == "C" & exp_num %in% c("M3","M4", "M5", "M8")| #Only use C mice from these experiments. Allocated groups to figures based on paper outline.
           group == "WM" & exp_num %in% c("M3","M4", "M5", "M8")|
           group == "WMC" & exp_num %in% c("M3","M4")|
           group == "WMR" & exp_num %in% c("M5","M6")) %>%
  mutate(group=factor(group, levels=c("C", "WM", "WMC", "WMR"))) 


#Post CDI PEG Subset
post_cdi_PEG_metadata <- metadata %>%
filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon")) %>% #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected
  filter(group == "C" & exp_num %in% c("M7","M9")| #Only use C mice from these experiments. Allocated groups to figures based on paper outline.
           group == "CWM" & exp_num %in% c("M6","M7", "M9")|
           group == "RM" & exp_num %in% c("M7","M9")|
           group == "FRM" & exp_num %in% c("M9")) %>%
  mutate(group=factor(group, levels=c("C", "CWM", "RM", "FRM")))  

#inner join to check for overlap
one_day_PEG %>% rbind(five_day_PEG) %>% inner_join(post_cdi_PEG)


#Split sequenced samples into 1 day PEG, 5 day PEG, and post-CDI subsets----------

#Full list of samples with sequence data that includes all timepoints
pcoa_data <- read_tsv("data/process/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>% slice_sample(n=100) %>%
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

