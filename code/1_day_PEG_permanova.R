source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Loads 16S output files for all sequenced samples

set.seed(19760620) #Same seed used for mothur analysis

#Remove variables not needed for this analysis
rm(agg_otu_data, diversity_data)

metadata <- metadata %>% 
  mutate(day = factor(day, levels = c(unique(as.factor(day)), "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0"))) #Transform day variable into factor variable

#PERMANOVA of 1-day PEG subset----
one_d <- read_dist("data/process/1_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
one_d_variables <- tibble(unique_label = attr(one_d, "Labels")) %>%
  left_join(metadata, by = "unique_label")
one_d_adonis <- adonis(one_d~(group/(exp_num*unique_cage_no*ext_plate*miseq_run))*day, strata = one_d_variables$unique_mouse_id, data = one_d_variables, permutations = 1000, parallel = 32)
one_d_adonis_table <- as_tibble(rownames_to_column(one_d_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/1_day_PEG_permanova.tsv")