source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Loads 16S output files for all sequenced samples

set.seed(19760620) #Same seed used for mothur analysis

#Remove variables not needed for this analysis
rm(agg_otu_data, diversity_data)

metadata <- metadata %>% 
  mutate(day = factor(day, levels = c(unique(as.factor(day)), "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0"))) #Transform day variable into factor variable

#PERMANOVA of 5-days PEG subset----
#Stools
five_d_stools <- read_dist("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
five_d_stools_variables <- tibble(unique_label = attr(five_d_stools, "Labels")) %>%
  left_join(metadata, by = "unique_label")
five_d_stools_adonis <- adonis(five_d_stools~(group/(exp_num*unique_cage_no*ext_plate*miseq_run))*day, strata = five_d_stools_variables$unique_mouse_id, data = five_d_stools_variables, permutations = 1000, parallel = 32)
five_d_stools_adonis_table <- as_tibble(rownames_to_column(five_d_stools_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/5_day_PEG_permanova_stools.tsv")
#Tissues
five_d_tissues <- read_dist("data/process/5_day_PEG/tissues/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
five_d_tissues_variables <- tibble(unique_label = attr(five_d_tissues, "Labels")) %>%
  left_join(metadata, by = "unique_label")
#PERMANOVA of all relevant variables (exclude Miseq run since they were the same for all samples in this subset)
five_d_tissues_adonis <- adonis(five_d_tissues~(group/(sample_type*exp_num*unique_cage_no*ext_plate))*day, strata = five_d_tissues_variables$unique_mouse_id, data = five_d_tissues_variables, permutations = 1000, parallel = 32)
five_d_tissues_adonis
#Write PERMANOVA results to tsv
five_d_tissues_adonis_table <- as_tibble(rownames_to_column(five_d_tissues_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/5_day_PEG_permanova_tissues.tsv")

#PERMANOVA of 1-day PEG subset----
one_d <- read_dist("data/process/1_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
one_d_variables <- tibble(unique_label = attr(one_d, "Labels")) %>%
  left_join(metadata, by = "unique_label")
one_d_adonis <- adonis(one_d~(group/(exp_num*unique_cage_no*ext_plate*miseq_run))*day, strata = one_d_variables$unique_mouse_id, data = one_d_variables, permutations = 1000, parallel = 32)
one_d_adonis_table <- as_tibble(rownames_to_column(one_d_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/1_day_PEG_permanova.tsv")

#PERMANOVA of post-CDI PEG subset----
#Stools
post_stools <- read_dist("data/process/post_CDI_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")
post_stools_variables <- tibble(unique_label = attr(post_stools, "Labels")) %>%
  left_join(metadata, by = "unique_label")
post_stools_adonis <- adonis(post_stools~(group/(exp_num*unique_cage_no*ext_plate*miseq_run))*day, strata = post_stools_variables$unique_mouse_id, data = post_stools_variables, permutations = 1000, parallel = 32)
post_stools_adonis_table <- as_tibble(rownames_to_column(post_stools_adonis$aov.tab, var = "effects")) %>% 
  write_tsv("data/process/post_CDI_PEG_permanova_stools.tsv")
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

#PERMANOVA of all samples together
# Read in Bray-Curtis distance matrix that represents all sequenced samples
all_dist <- read_dist("data/process/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata and prep_variables in order to test the impact of the relevant variables on community structure
all_variables <- tibble(unique_label = attr(all_dist, "Labels")) %>%
  left_join(metadata, by = "unique_label")

# See how plate number and library number contribute to community structure----
#Check subset of variables
#Test for sample_type, ext_plate, and any interactions between them.
sample_miseq_ext_plate_adonis <- adonis(all_dist~ sample_type*(miseq_run/ext_plate), data = all_variables, permutations = 1000, parallel = 32)
#Write PERMANOVA (adonis) results to .tsv
tibble(effects = c("sample_type", "miseq_run", "miseq_run:ext_plate", "sample_type:miseq_run:ext_plate", "Residuals"),
       r_sq = sample_miseq_ext_plate_adonis$aov.tab$R2[1:5],
       p = sample_miseq_ext_plate_adonis$aov.tab$Pr[1:5]) %>%
  write_tsv("exploratory/notebook/adonis_sample_miseq_plate.tsv")

#Plot PCoA data with different variables of interest assigned to the color aesthetic----
plot_pcoa_variable <- function(color_variable){
  all_pcoa_data %>%
    left_join(prep_variables, by = "unique_label") %>% #Add variables for ext plate number and miseq run
    ggplot(aes(x=axis1, y=axis2, color = {{ color_variable }}))+
    geom_point(size=2)+
    theme_classic()+
    theme(legend.position = "bottom")+
    xlim(-0.425, 0.65)+
    ylim(-0.525, 0.5)+
    labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
         y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))
}

#PCoA plot by sample_type
plot_pcoa_variable(sample_type)+
  ggsave("exploratory/notebook/pcoa_sample_type.pdf")
plot_pcoa_variable(ext_plate)+
  ggsave("exploratory/notebook/pcoa_ext_plate.pdf")
plot_pcoa_variable(miseq_run)+
  ggsave("exploratory/notebook/pcoa_miseq_run.pdf")

#PERMANOVA (adonis) of all relevant variables----
#Comment out- takes too much memory
#all_adonis <- adonis(all_dist~(group/(sample_type*exp_num*unique_cage_no*ext_plate*miseq_run))*day, strata = all_variables$unique_mouse_id, data = all_variables, permutations = 1000, parallel = 32)
#How to account for sample_type?
#Error: vector memory exhausted (limit reached?)
