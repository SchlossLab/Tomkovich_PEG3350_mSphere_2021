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
library(BiocManager)
#To install additional packages required for MetaLonDA with BioConductor
#BiocManager::install(c("metagenomeSeq", "DESeq2", "edgeR"))
#Will also ask you to update older versions of packages already installed. For now, select n for update none
library(MetaLonDA)

#Load read counts of 8 features from 100 samples
#Samples are from 2 groups, 5 subjects per group, and 10 time points per subject
data("metalonda_test_data")
#Create Group, Time & ID annotation vectors
n.group = 2
n.sample = 5
n.timepoints = 10
Group = factor(c(rep("A", n.sample*n.timepoints), rep("B", n.sample*n.timepoints)))
Time = rep(rep(1:n.timepoints, times = n.sample), 2)
ID = factor(rep(1:(2*n.sample), each = n.timepoints))

#Define the prediction timepoints
points = seq(1, 10, length.out = 100)
output.metalonda.f5 = metalonda(Count = metalonda_test_data[5,], Time = Time, Group = Group, ID = ID,
                                n.perm = 20, fit.method = "nbinomial", points = points, 
                                text = rownames(metalonda_test_data)[5], parall = FALSE, pvalue.threshold = 0.5,
                                adjust.method = "BH", time.unit = "hours", ylabel = "Read Counts", col = c("chartreuse", "blue4"))
#Note for test example 20 permutations were used, for real analysis. Recommend at least 1000 permutations

## Identify significant time intervals for all features: 
output.metalonda.all = metalondaAll(Count = metalonda_test_data, Time = Time, Group = Group,
                                    ID = ID, n.perm = 100, fit.method = "nbinomial", num.intervals = 100, 
                                    parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "hours", 
                                    norm.method = "none", prefix = "Test_metalondaALL", ylabel = "Read Counts",
                                    col = c("black", "green"))

source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files
#Remove unneeded files
rm(all_axis_labels, all_pcoa_data, diversity_data, pcoa_tissues_v_stool, pcoa_water, seq_prep_metadata,
   taxa_info, taxonomy)

#MetaLonDA can only compare 2 groups at a time
#Create dataframe of mice from WM and C groups and only stool samples
wm_v_c <- agg_otu_data %>% 
  filter(group %in% c("WM", "C") & sample_type == "stool" & day > 0) %>% #Only interested in modeling after perturbation
  select(unique_label, group, day, unique_mouse_id, otu, agg_rel_abund) %>% 
  mutate(otu = str_replace_all(otu, "/", "-")) #Remove / which cause issues with metalonda analysis

#remove otus that have all 0s across samples
wm_v_c_proc <- wm_v_c %>% 
  pivot_wider(id_cols = unique_label, names_from = otu, values_from = agg_rel_abund)  %>% 
  select(where(~ any(. !=0))) %>% #removes all columns with zero
  gather(-unique_label, key="otu", value="agg_rel_abund") %>% #Narrow dataframe
  pivot_wider(id_cols = otu, names_from = unique_label, values_from = agg_rel_abund) %>% #Widen with samples across top and otus as rows
  column_to_rownames(var = "otu") #Transform otu column into row names

wm_v_c_variables <- wm_v_c %>% 
  pivot_wider(id_cols = unique_label, names_from = otu, values_from = agg_rel_abund)  %>% 
  select(where(~ any(. !=0))) %>% #removes all columns with zero
  gather(-unique_label, key="otu", value="agg_rel_abund") %>% #Narrow dataframe
  pivot_wider(id_cols = otu, names_from = unique_label, values_from = agg_rel_abund) %>% #Widen with samples across top and otus as rows
  gather(-otu, key="unique_label", value="agg_rel_abund") %>% #Narrow dataframe 
  distinct(unique_label) %>% 
  left_join(metadata, by = "unique_label") %>% 
  mutate(day = as.integer(day))  #Day variable (transformed to integer to get rid of decimals on PCoA animation

group <- wm_v_c_variables$group
time <- wm_v_c_variables$day
id <- wm_v_c_variables$unique_mouse_id

#Figure out last timepoint in the dataset
max_day <- max(wm_v_c_variables$day)

#Group 1 color (determined alphabetically?)
group1_color <- "#238b45"
#Group 2 color (determined alphabetically?)
group2_color <- "#88419d"

#Define the prediction timepoints
points <- seq(1, max_day, length= max_day-1) 
output.metalonda.all <-  metalondaAll(Count = wm_v_c_proc, Time = time, Group = group,
                                    ID = id, n.perm = 1000, fit.method = "nbinomial", num.intervals = 46, 
                                    parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "days", 
                                    norm.method = "none", prefix = "Test_metalondaALL", ylabel = "Relative abundance",
                                    col = c(group1_color, group2_color))
#Time needed for 624 OTUs and 292 samples with 1000 permutations = 12 minutes per OTU
