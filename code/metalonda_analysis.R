#To install MetaLonDA
#install.packages("MetaLonDA")
#may need to change CRANmirror: chooseCRANmirror() Make a selection by typing in the number of desired option
#Also need to install the following packages that are dependencies:
#‘metagenomeSeq’, ‘DESeq2’, ‘edgeR’
#Not available for R version 4.0.2
#metagenomeSeq suggests installing Bioconductor packages with BiocManager
#To install BiocManager:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.12")
#Once BiocManager installed:
#library(BiocManager) #Only need to install packages
#To install additional packages required for MetaLonDA with BioConductor
#BiocManager::install(c("metagenomeSeq", "DESeq2", "edgeR"))
#Will also ask you to update older versions of packages already installed. For now, select n for update none
library(MetaLonDA)

#Load read counts of 8 features from 100 samples
#Samples are from 2 groups, 5 subjects per group, and 10 time points per subject
#data("metalonda_test_data")
#Create Group, Time & ID annotation vectors
# n.group = 2
# n.sample = 5
# n.timepoints = 10
# Group = factor(c(rep("A", n.sample*n.timepoints), rep("B", n.sample*n.timepoints)))
# Time = rep(rep(1:n.timepoints, times = n.sample), 2)
# ID = factor(rep(1:(2*n.sample), each = n.timepoints))

#Define the prediction timepoints
# points = seq(1, 10, length.out = 100)
# output.metalonda.f5 = metalonda(Count = metalonda_test_data[5,], Time = Time, Group = Group, ID = ID,
#                                 n.perm = 20, fit.method = "nbinomial", points = points,
#                                 text = rownames(metalonda_test_data)[5], parall = FALSE, pvalue.threshold = 0.5,
#                                 adjust.method = "BH", time.unit = "hours", ylabel = "Read Counts", col = c("chartreuse", "blue4"))
# #Note for test example 20 permutations were used, for real analysis. Recommend at least 1000 permutations
#
# ## Identify significant time intervals for all features:
# output.metalonda.all = metalondaAll(Count = metalonda_test_data, Time = Time, Group = Group,
#                                     ID = ID, n.perm = 100, fit.method = "nbinomial", num.intervals = 100,
#                                     parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "hours",
#                                     norm.method = "none", prefix = "Test_metalondaALL", ylabel = "Read Counts",
#                                     col = c("black", "green"))

source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files
#Remove unneeded files
rm(all_axis_labels, all_pcoa_data, diversity_data, pcoa_tissues_v_stool, pcoa_water, seq_prep_metadata,
   taxa_info, taxonomy)

#Make a function to use Metalonda for the PEG3350 project
#MetaLonDA can only compare 2 groups at a time (list groups alpabetically)
#Function is set up to look at timepoints after day 0 of the experiment (model fitting method does not work well for large perturbations)
#Arguments:
#group1 <- name of 1st group in dataset to compare
#group2 <- name of 1st group in dataset to compare
#specify_sample_type <- type of samples to analyze, if more than one enclose in c()
#group1_color <- color to use to plot group 1
#group2_color <- color to use to plot group 2
#no_permutations <- number of permutations to run
run_metalonda <- function(group1, group2, specify_sample_type, group1_color, group2_color, permutations){
  #Create dataframe of mice from the 2 groups and only sample_type specified
  sample_list <- agg_otu_data %>%
    filter(group %in% c(group1, group2) & sample_type %in% "stool" & day > 0) %>% #Only interested in modeling after perturbation
    select(unique_label, group, day, unique_mouse_id, otu, agg_rel_abund) %>%
    mutate(otu = str_replace_all(otu, "/", "-")) %>% #Remove / which cause issues with metalonda analysis
    pivot_wider(id_cols = unique_label, names_from = otu, values_from = agg_rel_abund)  %>%
    select(where(~ any(. !=0)))  #removes all columns with zero

  #remove otus that have all 0s across samples
  processed_samples <- sample_list %>%
    gather(-unique_label, key="otu", value="agg_rel_abund") %>% #Narrow dataframe
    pivot_wider(id_cols = otu, names_from = unique_label, values_from = agg_rel_abund) %>% #Widen with samples across top and otus as rows
    column_to_rownames(var = "otu") #Transform otu column into row names

  sample_variables <- sample_list %>%
    gather(-unique_label, key="otu", value="agg_rel_abund") %>% #Narrow dataframe
    pivot_wider(id_cols = otu, names_from = unique_label, values_from = agg_rel_abund) %>% #Widen with samples across top and otus as rows
    gather(-otu, key="unique_label", value="agg_rel_abund") %>% #Narrow dataframe
    distinct(unique_label) %>%
    left_join(metadata, by = "unique_label") %>%
    mutate(day = as.integer(day))  #Day variable (transformed to integer to get rid of decimals on PCoA animation

  group <- sample_variables$group
  time <- sample_variables$day
  id <- sample_variables$unique_mouse_id

  #Figure out last timepoint in the dataset
  max_day <- max(sample_variables$day)

  #Define the prediction timepoints
  points <- seq(1, max_day, length= max_day-1)
  output.metalonda.all <-  metalondaAll(Count = processed_samples, Time = time, Group = group,
                                        ID = id, n.perm = permutations, fit.method = "nbinomial", num.intervals = 46,
                                        parall = TRUE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "days",
                                        norm.method = "none", prefix = "Test_metalondaALL", ylabel = "Relative abundance",
                                        col = c(group1_color, group2_color))

}
#parall argument = TRUE in metalondaAll will detect max cores available and use 1 - the max cores.


#Time needed per OTU = 14 minutes

#Run metalonda for C & WM group stool samples and 1000 permutations
#c_vs_wm <- run_metalonda("C", "WM", "stool", "#238b45", "#88419d", 1000)
#Move output of metalondaAll to new folder before running again
# mkdir Test_metalondaCvsWM 
#mv Test_metalondaALL/* Test_metalondaCvsWM/

#Run metalonda for M1 & WM group stool samples and 1000 permutations
#Create sequential color scheme to compare different PEG groups since we've reused colors
#Help picking color scheme: https://hihayk.github.io/scale/#4/6/50/80/-51/67/20/14/88419d/136/65/157/white
m1_vs_wm <- run_metalonda("M1", "WM", "stool", "#BDC3E1", "#88419d", 1000)

#978AC7 post-CDI PEG
