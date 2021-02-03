source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Loads 16S output files for all sequenced samples

set.seed(19760620) #Same seed used for mothur analysis

#Remove variables not needed for this analysis
rm(agg_otu_data, diversity_data)

metadata <- metadata %>% 
  mutate(day = factor(day, levels = c(unique(as.factor(day)), "not_applicable", "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0"))) #Transform day variable into factor variable

#PERMANOVA of post-CDI PEG subset----
#Tissues
#Tissues
post_tissues <- read_dist("data/process/post_CDI_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
post_tissues_variables <- tibble(unique_label = attr(post_tissues, "Labels")) %>%
  left_join(metadata, by = "unique_label")
#PERMANOVA of all relevant variables (exclude Miseq run since they were the same for all samples in this subset)
post_tissues_adonis <- adonis(post_tissues~sample_type*unique_cage_no, strata = post_tissues_variables$unique_mouse_id, data = post_tissues_variables, permutations = 1000, parallel = 32)
post_tissues_adonis
#Write PERMANOVA results to tsv
post_tissues_adonis_table <- as_tibble(rownames_to_column(post_tissues_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_tissues.tsv")

#Stools
post_stools <- read_dist("data/process/post_CDI_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
post_stools_variables <- tibble(unique_label = attr(post_stools, "Labels")) %>%
  left_join(metadata, by = "unique_label")

#Look at variables of interest individually----
#group
group_adonis <- adonis(post_stools~group, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
group_adonis_table <- as_tibble(rownames_to_column(group_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_group.tsv")
#day
day_adonis <- adonis(post_stools~day, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
day_adonis_table <- as_tibble(rownames_to_column(day_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_day.tsv")
#exp_num
exp_num_adonis <- adonis(post_stools~exp_num, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
exp_num_adonis_table <- as_tibble(rownames_to_column(exp_num_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_exp_num.tsv")
#unique_cage_no
cage_adonis <- adonis(post_stools~unique_cage_no, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
cage_adonis_table <- as_tibble(rownames_to_column(cage_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_cage.tsv")
#ext_plate
ext_plate_adonis <- adonis(post_stools~ext_plate, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
ext_plate_adonis_table <- as_tibble(rownames_to_column(ext_plate_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_ext_plate.tsv")
#miseq_run
miseq_run_adonis <- adonis(post_stools~miseq_run, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
miseq_run_adonis_table <- as_tibble(rownames_to_column(miseq_run_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_miseq_run.tsv")
#g_d = group * day
g_d_adonis <- adonis(post_stools~group*day, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
g_d_adonis_table <- as_tibble(rownames_to_column(g_d_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_group_day.tsv")
#g_c = group * unique_cage_no
g_c_adonis <- adonis(post_stools~group * unique_cage_no, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
g_c_adonis_table <- as_tibble(rownames_to_column(g_c_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_group_cage.tsv")
#g_e = group * exp_num
g_e_adonis <- adonis(post_stools~group * exp_num, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
g_e_adonis_table <- as_tibble(rownames_to_column(g_e_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_group_exp_num.tsv")
#g_ex = group * ext_plate
g_ex_adonis <- adonis(post_stools~group * ext_plate, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
g_ex_adonis_table <- as_tibble(rownames_to_column(g_ex_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_group_ext_plate.tsv")
#g_m = group * miseq_run
g_m_adonis <- adonis(post_stools~group * miseq_run, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 10)
g_m_adonis_table <- as_tibble(rownames_to_column(g_m_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools_group_miseq_run.tsv")
#Look at majority of variables without nesting---
post_stools_adonis <- adonis(post_stools~(group/(unique_cage_no*ext_plate))*day, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 10, parallel = 10)
post_stools_adonis_table <- as_tibble(rownames_to_column(post_stools_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools.tsv")


#Initial PERMANOVA design with all variables of interest was taking too long to run
#~7 days for 10 iterations did not finish running. So changed the design (see above) to examine variables of interest individually
#post_stools_adonis <- adonis(post_stools~(group/(sample_type*exp_num*unique_cage_no*ext_plate*miseq_run))*day, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 10, parallel = 10)
#post_stools_adonis_table <- as_tibble(rownames_to_column(post_stools_adonis$aov.tab, var = "effects")) %>% 
#  write_tsv("data/process/post_CDI_PEG_permanova_stools.tsv")
