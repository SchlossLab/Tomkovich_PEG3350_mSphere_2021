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
five_d_tissues_adonis <- adonis(five_d_tissues~(group/(sample_type*exp_num*unique_cage_no*ext_plate))*day, strata = five_d_tissues_variables$unique_mouse_id, data = five_d_tissues_variables, permutations = 1000, parallel = 10)
five_d_tissues_adonis
#Write PERMANOVA results to tsv
five_d_tissues_adonis_table <- as_tibble(rownames_to_column(five_d_tissues_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_tissues.tsv")

#Stools
five_d_stools <- read_dist("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
five_d_stools_variables <- tibble(unique_label = attr(five_d_stools, "Labels")) %>%
  left_join(metadata, by = "unique_label")

#Look at variables of interest individually
#group
group_stools_adonis <- adonis(five_d_stools~group, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
group_stools_adonis_table <- as_tibble(rownames_to_column(group_stools_adonis$aov.tab, var = "effects")) %>%
write_tsv("data/process/5_days_PEG_permanova_stools_group.tsv")
#day
day_stools_adonis <- adonis(five_d_stools~day, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
day_stools_adonis_table <- as_tibble(rownames_to_column(day_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_day.tsv")
#exp_num
exp_num_stools_adonis <- adonis(five_d_stools~exp_num, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
exp_num_stools_adonis_table <- as_tibble(rownames_to_column(exp_num_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_exp_num.tsv")
#unique_cage_no
cage_stools_adonis <- adonis(five_d_stools~unique_cage_no, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
cage_stools_adonis_table <- as_tibble(rownames_to_column(cage_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_cage.tsv")
#ext_plate
ext_plate_stools_adonis <- adonis(five_d_stools~ext_plate, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
ext_plate_stools_adonis_table <- as_tibble(rownames_to_column(ext_plate_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_ext_plate.tsv")
#miseq_run
miseq_run_stools_adonis <- adonis(five_d_stools~miseq_run, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
miseq_run_stools_adonis_table <- as_tibble(rownames_to_column(miseq_run_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_miseq_run.tsv")
#g_d = group * day
g_d_stools_adonis <- adonis(five_d_stools~group*day, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
g_d_stools_adonis_table <- as_tibble(rownames_to_column(g_d_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_group_day.tsv")
#g_c = group * unique_cage_no
g_c_stools_adonis <- adonis(five_d_stools~group*unique_cage_no, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
g_c_stools_adonis_table <- as_tibble(rownames_to_column(g_c_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_group_cage.tsv")
#g_e = group * exp_num
g_e_stools_adonis <- adonis(five_d_stools~group*exp_num, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
g_e_stools_adonis_table <- as_tibble(rownames_to_column(g_e_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_group_exp_num.tsv")
#g_ex = group * ext_plate
g_ex_stools_adonis <- adonis(five_d_stools~group*ext_plate, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
g_ex_stools_adonis_table <- as_tibble(rownames_to_column(g_ex_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_group_ext_plate.tsv")
#g_m = group * miseq_run
g_m_stools_adonis <- adonis(five_d_stools~group*miseq_run, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 10)
g_m_stools_adonis_table <- as_tibble(rownames_to_column(g_m_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools_group_miseq_run.tsv")
#Look at majority of variables without nesting---
five_d_stools_adonis <- adonis(five_d_stools~day*group*unique_cage_no*ext_plate*miseq_run*exp_num, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 10, parallel = 10)
five_d_stools_adonis_table <- as_tibble(rownames_to_column(five_d_stools_adonis$aov.tab, var = "effects")) %>%
  write_tsv("data/process/5_days_PEG_permanova_stools.tsv")

#Initial PERMANOVA design with all variables of interest was taking too long to run
#~7 days for 10 iterations did not finish running. So changed the design (see above) to examine variables of interest individually
#five_d_stools_adonis <- adonis(five_d_stools~(group/(exp_num*unique_cage_no*ext_plate*miseq_run))*day, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 10, parallel = 10)
#five_d_stools_adonis_table <- as_tibble(rownames_to_column(five_d_stools_adonis$aov.tab, var = "effects")) %>%
#  write_tsv("data/process/5_days_PEG_permanova_stools.tsv")
