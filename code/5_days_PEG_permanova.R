source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Loads 16S output files for all sequenced samples

set.seed(19760620) #Same seed used for mothur analysis

#Remove variables not needed for this analysis
rm(agg_otu_data, diversity_data)

metadata <- metadata %>% 
  mutate(day = factor(day, levels = c(unique(as.factor(day)), "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0"))) #Transform day variable into factor variable

#PERMANOVA of 5-days PEG subset----
#Tissues
five_d_tissues <- read_dist("data/process/5_day_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
five_d_tissues_variables <- tibble(unique_label = attr(five_d_tissues, "Labels")) %>%
  left_join(metadata, by = "unique_label")
#PERMANOVA of all relevant variables (exclude Miseq run since they were the same for all samples in this subset)
five_d_tissues_adonis <- adonis(five_d_tissues~(group/(sample_type*exp_num*unique_cage_no*ext_plate))*day, strata = five_d_tissues_variables$unique_mouse_id, data = five_d_tissues_variables, permutations = 1000, parallel = 32)
five_d_tissues_adonis
#Write PERMANOVA results to tsv
five_d_tissues_adonis_table <- as_tibble(rownames_to_column(five_d_tissues_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/5_days_PEG_permanova_tissues.tsv")

#Stools
five_d_stools <- read_dist("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
five_d_stools_variables <- tibble(unique_label = attr(five_d_stools, "Labels")) %>%
  left_join(metadata, by = "unique_label")
five_d_stools_adonis <- adonis(five_d_stools~(group/(exp_num*unique_cage_no*ext_plate*miseq_run))*day, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 32)
five_d_stools_adonis_table <- as_tibble(rownames_to_column(five_d_stools_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/5_days_PEG_permanova_stools.tsv")

