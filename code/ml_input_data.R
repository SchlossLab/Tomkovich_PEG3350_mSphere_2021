source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Goal: use Machine learning models to predict mice that clear within 10 days or remain colonized for longer than 10 days

#For this analysis, we will only examine the groups of mice that were challenged with C. diff
cdiff_groups <- c("C", "WM", "WMC", "WMR",
                "M1", "1RM1",
                "CWM", "FRM", "RM")

#Figure out the timepoints that have the largest number of samples with sequencing data
#Select samples from only mice that were challenged with C. difficile
cdiff_diversity <- diversity_data %>% 
  filter(group %in% cdiff_groups)
#Select stool samples only
cdiff_diversity_stools <- subset_stool(cdiff_diversity)

#Examine how many stool samples we have per timepoint
num_days <- cdiff_diversity_stools %>% 
  group_by(day) %>%
  count() %>% 
  arrange(desc(n))

#Choose to examine input data from D5 & D1 since these days have the largest number of samples associated with them

#Select metadata rows where we no d10 C. difficile colonization status
metadata <- metadata %>% 
  filter(!clearance_status_d10 == "no_data")

#Read in OTU data----
otu_data <- read_tsv("data/process/peg3350.opti_mcc.0.03.subsample.shared", col_types=cols(Group=col_character())) %>%
  select(-label, -numOtus) %>%
  rename(unique_label = Group) %>% #group is the same as unique_label in the metadata data frame
  gather(-unique_label, key="otu", value="count") %>%
  mutate(rel_abund=count/1000) %>%  #Use 1000, because this is the subsampling parameter chosen.
  filter(!otu == "Otu00012") %>% #Remove C. difficile OTU
  pivot_wider(id_cols = unique_label, names_from = otu, values_from = rel_abund)%>% #Transform dataframe so that each OTU is a different column
  left_join(select(metadata, clearance_status_d10, unique_label, day), by = "unique_label") %>% 
  rename(outcome = clearance_status_d10) %>%  #Rename group to outcome, this is what we will classify based on OTU relative abundances
  relocate(outcome, before = unique_label) 

#Function to subset OTU data to a specific timepoint  
subset_ml_input <- function(timepoint){
  otu_data %>% 
    filter(day == timepoint) %>% 
    select(-unique_label, -day) %>% #drop columns we no longer need
    write_csv(path = paste0("data/process/ml_input_d", timepoint, ".csv"))
}

#Create 5dpi ml input data to predict d10 clearance status based on bacterial communities-----
d5_input <- subset_ml_input(5)
#77 mice
#Examine the balance of cleared vs colonized mice
d5_input %>% count(outcome)
#37 cleared, 40 mice colonized
#Create 1dpi ml input data to predict d10 clearance status based on bacterial communities-----
d1_input <- subset_ml_input(1)
#70 mice
d1_input %>% count(outcome)
#38 mice cleared, 32 mice colonized

#Consider making additional input data with WMR group, seems more unique in terms of being able to detect
#C. difficile relative abundance in the 16S data compared to the other PEG groups