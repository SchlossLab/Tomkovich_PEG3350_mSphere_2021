source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Goal: use Machine learning models to predict mice that clear within 10 days or remain colonized for longer than 10 days

#For this analysis, we will only examine the groups of mice that were challenged with C. diff
cdiff_groups <- c("C", "WM", "WMC",
                "M1", "1RM1",
                "CWM", "FRM", "RM")
#Remove WMR group because appearance in C. diff was delayed relative to the other PEG-treated mice

#Figure out the timepoints that have the largest number of samples with sequencing data
#Select samples from only mice that were challenged with C. difficile
#Read in alpha diversity metrics for all samples----
diversity_data <- read_tsv("data/process/peg3350.opti_mcc.groups.ave-std.summary") %>%
  filter(method == "ave") %>%
  select(group, sobs, shannon, invsimpson, coverage) %>%
  rename(unique_label = group) %>% #group is the same as unique_label in the metadata data frame
  left_join(metadata, by = "unique_label") #Match only the samples we have sequence data for

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
  filter(!clearance_status_d10 == "no_data") %>% 
  filter(!group == "WMR") #Remove WMR group

#Read in OTU data----
otu_data <- read_tsv("data/process/peg3350.opti_mcc.0.03.subsample.shared", col_types=cols(Group=col_character())) %>%
  select(-label, -numOtus) %>%
  rename(unique_label = Group) %>% #group is the same as unique_label in the metadata data frame
  gather(-unique_label, key="otu", value="count") %>%
  mutate(rel_abund=count/1000) %>%  #Use 1000, because this is the subsampling parameter chosen.
  filter(!otu == "Otu0012") %>% #Remove C. difficile OTU
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
#69 mice
#Examine the balance of cleared vs colonized mice
d5_input %>% count(outcome)
#36 cleared, 33 mice colonized
#Create 1dpi ml input data to predict d10 clearance status based on bacterial communities-----
d1_input <- subset_ml_input(1)
#61 mice
d1_input %>% count(outcome)
#36 mice cleared, 25 mice colonized

#Excluded WMR group because these mice seemed more unique in terms of being able to detect C. diff
#C. difficile relative abundance in the 16S data compared to the other PEG groups