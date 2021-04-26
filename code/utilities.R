#Load required libraries
library(tidyverse)
library(readxl)
library(writexl)
library(cowplot)
library(broom)
library(knitr)
library(ggpubr)
library(gganimate)
library(glue)
library(ggtext)
library(vegan)
library(parallel)
library(magick)

#Read in metadata----
metadata <- read_excel("data/process/metadata.xlsx", col_types = c("text", "numeric", "text", "text", "numeric", "text", "numeric", "text", "text", "text", "text", "numeric", "numeric")) #specify column types

#Check for duplicated unique_labels
duplicated <- metadata %>%
  filter(duplicated(unique_label)) #0 duplicates

#Read in metadata from library preparation for samples that underwent 16S rRNA gene sequencing----
seq_prep_metadata <- read_tsv("data/process/16Sprep_PEG3350_metadata", col_types = "ifffdcfcfffDDDfDf") %>% #specify the col_types]
  mutate(unique_label = replace(unique_label, unique_label == "FMT_from_Motility_9", "FMTMotility9"), #rename FMT & PBS gavage samples by removing the underscores to match fastq file names 
         unique_label = replace(unique_label, unique_label == "PBS_Gavage_from_D4_Motility_9", "PBSD4Motility9"))

#Identify samples that were sequenced twice in metadata
duplicated <- seq_prep_metadata %>%
  filter(duplicated(unique_label)) %>% #5 samples sequenced twice
  pull(unique_label)

#Metadata for samples that were sequenced twice:
duplicated_seq_samples <- seq_prep_metadata %>%
  filter(unique_label %in% duplicated)

#The notes column has important notes from the 16S library preparation that should be examined before 16S rRNA sequencing analysis:
prep_notes <- unique(seq_prep_metadata$notes)
#13 different types of notes, rest are NA (blank)
#notes that indicate contaminated samples that should be dropped from the dataset:
"NOTE from Lucas: columns 9 and 10 contain material from 9-12 on the bead plate because of epMotion mistake."
"NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically."
"Lucas' note had a C or 6, left sample label as is, NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically."
"Lucas' note had a WM 18 listed, but there is no mouse with this group and mouse number label from motility, kept as is, NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically."
"Lucas' note had day 7 but likely a legability issue since the samples are organized chronologically in the boxes., NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically."
"Leftover stool from when this sample was arrayed for sequencing in plate_2, NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically."
#List of values in notes that indicate samples were contaminated:
contaminated_notes <- c("NOTE from Lucas: columns 9 and 10 contain material from 9-12 on the bead plate because of epMotion mistake.",
                        "NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically.",
                        "Lucas' note had a C or 6, left sample label as is, NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically.",
                        "Lucas' note had a WM 18 listed, but there is no mouse with this group and mouse number label from motility, kept as is, NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically.",
                        "Lucas' note had day 7 but likely a legability issue since the samples are organized chronologically in the boxes., NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically.",
                        "Leftover stool from when this sample was arrayed for sequencing in plate_2, NOTE from Lucas: Columns 3,4,5 contain contamination from plate_7, column 1 due to epMotion error. A conservative estimate because only column 4 got the excess liquid from another column on plate 7, but since it was overflowing I couldn't tell if the contents ran into the columns next to it (3 and 5) or just stayed in a puddle in the grooves in the top of the plate. But better safe than sorry I figured. Check if possiblee to examine bioinformatically."
)

#Create list of unique_label for the contaminated files to remove
contaminated_samples <- seq_prep_metadata %>%
  filter(notes %in% contaminated_notes) %>%
  mutate(grep_prefix = "'^", # ^ indicates the start of a string
         grep_suffix = "_.*'") %>%  #Add columns to construct regexp to remove these sequences from the project data/raw folder to exclude contaminated samples from 16S gene rRNA sequencing analysis. . in _.* is to avoid matching to files that don't have an underscore after pattern is matched
  unite(col = "grep_name", grep_prefix, unique_label, grep_suffix, sep = "", remove = TRUE) %>% #merge the 3 columns into 1 column to create the bash commands
  pull(grep_name) %>%
  noquote() # remove quotes around bash commands
paste(contaminated_samples, sep =" \n", collapse = " ") # print all rows
#Use this list with a for loop from the command line to remove these sequences from data raw
#See code/copy_fastqs_to_data

#Check samples with these notes during analysis:
"From Lucas: plates 5 and 6 were incorrectly labeled at 3 and 4 (I checked the dates on them to make sure they were the right plates) but because of the mix up it is possible that 5 and 6 may be flipped (small chance). So just double check that the results look ok,"
"Whole cecum frozen, Lucas took a snip for sequencing with scalpel"
"Note that sample yielded no DNA from Lucas"

#notes that need no further action (already corrected in file):
"5/17/19 was not labeleled on the side of the tube"
"Katie's note lists D19 for these samples instead of D20, which was already sequenced as part of plates1_2 library"
"swapped with A6 when samples were arrayed into extraction plates"
"swapped with B7 when samples were arrayed into extraction plates"
"Realized from Katie's note there was tube_label date typo that has now been corrected in make_metadata_file, D5 was initially entered as 3/21/19 instead of 3/31/19"

#Read in peg3350.files and cross check with seq_prep_metadata making sure there is no typos in NIAIDS, and all IDs match)
#NOTE: Samples from plates 14-17 have not yet been sequenced so are not in peg3350.files

peg3350.files <- read_tsv("data/process/peg3350.files", col_names=c("unique_label", "read_1", "read_2")) #no columns in .files format

peg3350.files_unique_label <- peg3350.files %>% select(unique_label) #Read only the unique label

#Check which samples in peg3350 files are not in seq_prep_metadata
seq_files_missing_from_metadata <- anti_join(peg3350.files, seq_prep_metadata) %>% 
  pull(unique_label)
#[1] "mock10"  "mock11"  "mock12"  "mock13"  "mock14"  "mock15"  "mock16"  "mock17" 
#[9] "mock1"   "mock2"   "mock3"   "mock4"   "mock5"   "mock6"   "mock7"   "mock8"  
#[17] "mock9"   "water10" "water11" "water12" "water13" "water14" "water15" "water16"
#[25] "water17" "water1"  "water2"  "water3"  "water4"  "water5"  "water6"  "water7" 
#[33] "water8"  "water9"
## All these control samples are in peg3350 files and not in seq_prep_metadata

#Removed the second set of sequences for all 5 duplicate's second sample from plate 8, 11, and 12 since all of these samples underwent an additional freeze thaw cycle
duplicates_to_drop <- seq_prep_metadata %>% 
  filter(unique_label %in% duplicated) %>% 
  filter(miseq_run != "motility_plates_1-2") #Drop all duplicate samples that were sequenced after the library that containted plates1-2

#Drop the 5 duplicates from seq_prep_metadata (these were already removed from data/raw, see code/copy_fastqs_to_data)
seq_prep_metadata <- seq_prep_metadata %>% 
  anti_join(duplicates_to_drop) %>% 
  select(unique_label, ext_plate, miseq_run) #Just select unique_label and the plate # and miseq run variables to test with adonis/PERMANOVA

#Unique_mice
metadata %>% filter(!duplicated(unique_mouse_id)) %>% #Filter to just unique mouse ids
  tally()
#150 mice total

#Make a column for cfu_d8 values
cfu_d8 <- metadata %>% 
  filter(day == "8") %>% 
  filter(!is.na(avg_cfu)) %>% #Remove NA values
  mutate(cfu_d8 = avg_cfu) %>% 
  select(unique_mouse_id, cfu_d8) #Will join to metadata by mouse id instead of unique_label so we can know each individual mouse's day 10 colonization status

#Make a column for cfu_d10 values
cfu_d10 <- metadata %>% 
  filter(day == "10") %>% 
  filter(!is.na(avg_cfu)) %>% #Remove NA values
  mutate(cfu_d10 = avg_cfu) %>% 
  select(unique_mouse_id, cfu_d10) #Will join to metadata by mouse id instead of unique_label so we can know each individual mouse's day 10 colonization status

#Join cfu_d8 & cfu_d10 columns to metadata
metadata <- metadata %>% 
  left_join(cfu_d8, by = "unique_mouse_id") %>% 
  left_join(cfu_d10, by = "unique_mouse_id")

#Join metadata to 16S seq prep metadata----
metadata <- seq_prep_metadata %>% 
  full_join(metadata, by = "unique_label") %>% 
  mutate(sample_type = case_when(str_detect(unique_label, "water") ~ "water",
                                 str_detect(unique_label, "FMT") ~ "FMT",
                                 TRUE ~ sample_type),
         group = case_when(str_detect(unique_label, "water") ~ "water",
                           str_detect(unique_label, "FMT") ~ "FMT",
                           TRUE ~ group)) %>% 
  #Replace NAs for any of the relevant variables (primarily just for the water and FMT samples)
  mutate(unique_mouse_id = replace_na(unique_mouse_id, "not_applicable")) %>% 
  mutate(unique_cage_no = replace_na(unique_cage_no, "not_applicable")) %>% 
  mutate(group = replace_na(group, "not_applicable")) %>% 
  mutate(exp_num = replace_na(exp_num, "not_applicable")) %>% 
  mutate(day = replace_na(day, "not_applicable")) %>% 
  #Make a column for subset designation (5-day, 1-day, post-CDI, or C for clindamycin mice)
  mutate(subset = case_when(group %in% c("WM", "WMR", "WMC", "WMN") ~ "5-day",
                            group %in% c("M1", "1RM1") ~ "1-day",
                            group %in% c("CWM", "FRM", "RM") ~ "post-CDI",
                            group %in% c("C", "CN") ~ "clind.",
                            TRUE ~ "not_applicable")) %>% 
  mutate(group = factor(group, levels = unique(as.factor(group))), #Transform group variable into factor variable
         unique_cage_no = factor(unique_cage_no, levels = unique(as.factor(unique_cage_no))), #Transform unique_cage_no variable into factor variable, add PT level to indicate pre-treatment samples
         exp_num = factor(exp_num, levels = unique(as.factor(exp_num))), #Transform exp_num variable into factor variable
         sample_type = factor(sample_type, levels = unique(as.factor(sample_type))), #Transform sample_type variable into factor variable
         ext_plate = factor(ext_plate, levels = unique(as.factor(ext_plate))), #Transform ext_plate variable into factor variable
         miseq_run = factor(miseq_run, levels = unique(as.factor(miseq_run))), #Transform miseq_run variable into factor variable
         unique_mouse_id = factor(unique_mouse_id, levels = unique(as.factor(unique_mouse_id)))) %>% #Transform unique_mouse_id into factor variable
  #Make a column for C. difficile clearance status at Day 10
  mutate(clearance_status_d10 = case_when(cfu_d10 > 0 ~ "colonized",
                                          cfu_d8 > 0 & cfu_d10 == 0 ~ "cleared",
                                          cfu_d8 == 0 & cfu_d10 > 0 ~ "colonized", #Some mice from WMR group did not show up as colonized until later timepoints, 10_M6 WMR mouse only showed up as colonized on d30?
                                          cfu_d8 == 0 ~ "cleared", #some of the mice from the 1-day PEG subset only have CFU quantified through d8
                                          cfu_d10 == 0 ~ "cleared",
                                          TRUE ~ "no_data")) %>% 
  #Make a column to denote whether mice were challenged with C. difficile (mock vs C. difficile challenged mice)
  mutate(infected = case_when(grepl("N", group) ~ "no", #Make a new column based on whether mice were challenged with C. difficile
                       TRUE ~ "yes")) #An N in the Group name indicates they were not challenged with C. difficile, all the other groups were challenged

#Check numbers of each group that have cleared or remain colonized with C. diff by d10
clear_v_col <- metadata %>% 
  group_by(group, clearance_status_d10) %>% 
  tally()
#A couple of mice from RM, FRM, CWM, and WMR groups were clear by 10.
#Some mice from groups of interest have no_data entries because they were only followed through d4 or d6 or died early

#Functions----

#Functions to define the 3 main subsets of mice used throughout the paper----
#Function to create 1 day PEG Subset of a given dataframe (df) :
one_day_PEG_subset <- function(df){
  df %>% 
    filter(group == "C" & exp_num %in% c("M6")| #Only use C mice from this experiments. Allocated groups to figures based on paper outline.
             group == "1RM1" & exp_num %in% c("M6R")| #Had to differentiate experiment 6 from 6R in the metadata to create unique_mouse_id that wouldn't overlap for the M1 & 1RM1 mice that are both labeled with mouse_ids that are #s1-6
             group == "M1" & exp_num %in% c("M6"))
} 

#5 days PEG Subset
five_day_PEG_subset <- function(df){
  df %>% 
    filter(group == "C" & exp_num %in% c("M3","M4", "M5", "M8")| #Only use C mice from these experiments. Allocated groups to figures based on paper outline.
             group == "WM" & exp_num %in% c("M3","M4", "M5", "M8")|
             group == "WMC" & exp_num %in% c("M3","M4")|
             group == "WMR" & exp_num %in% c("M5","M6")) 
}    

#Post CDI PEG Subset
post_cdi_PEG_subset <- function(df){
  df %>% 
    filter(group == "C" & exp_num %in% c("M7","M9")| #Only use C mice from these experiments. Allocated groups to figures based on paper outline.
             group == "CWM" & exp_num %in% c("M6","M7", "M9")|
             group == "RM" & exp_num %in% c("M7","M9")|
             group == "FRM" & exp_num %in% c("M9"))
} 

#Functions to split up the different sample type subsets within the 5-day PEG and Post-CDI PEG subsets----
#Function to filter dataframe (df) to examine stool samples
subset_stool <- function(df){
  df %>% filter(sample_type == "stool")
}
#Function to filter dataframe (df) to examine tissue samples
subset_tissue <- function(df){
  df %>% filter(!sample_type == "stool")
}

#Function to filter dataframe (df) to examine cecum tissue samples
subset_cecum <- function(df){
  df %>% filter(sample_type == "cecum")
}

#Function to filter dataframe (df) to examine proximal_colon tissue samples
subset_proximal_colon <- function(df){
  df %>% filter(sample_type == "proximal_colon")
}

#Function to filter dataframe (df) to examine distal_colon tissue samples
subset_distal_colon <- function(df){
  df %>% filter(sample_type == "distal_colon")
}

#Function to add the mock challenged mice (CN and WMN groups) to a dataframe (diversity_data or agg_otu_data)----
add_mocks <- function(subset_df, original_df){
  subset_df %>% 
    add_row(original_df %>% filter(group %in% c("CN", "WMN"))) #Add the groups of mock challenged mice
}

#Function to figure out how many samples we have per group per day for each subset (subset = dataframe of a subset of samples)----
count_subset <- function(subset){
  subset %>% 
    group_by(group) %>%
    count(day) %>% 
    arrange(day)
}

#Function to have y-axis in scientific notation----
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

#Functions used in statistical analysis----
#Function to calculate the median cfu values per group from a dataframe (x)
get_cfu_median <- function(x){
  x %>%
    group_by(group) %>%
    summarize(median=median(avg_cfu)) %>%
    spread(key=group, value=median)
}

#Function to calculate the median weight_change values per group from a dataframe (x)
get_weight_median <- function(x){
  x %>%
    group_by(group) %>%
    summarize(median=median(weight_change)) %>%
    spread(key=group, value=median)
}

#Function to calculate the median histology summary score values per group from a dataframe (x)
get_summary_score_median <- function(x){
  x %>%
    group_by(Group) %>%
    summarize(median=median(summary_score)) %>%
    spread(key=Group, value=median)
}

#Function to tidy pairwise comparisons to use for adding stats to plots
tidy_pairwise <- function(spread_pairwise){
  spread_pairwise %>%
    pivot_longer(-day, names_to = "compare", values_to = "p.adj") %>%
    separate(col = compare, c("group1", "group2"), sep = "-", remove = TRUE)
}
#Function to tidy pairwise comparisons to use for adding stats to otu plots----
tidy_pairwise_otu <- function(spread_pairwise){
  spread_pairwise %>%
    pivot_longer(-otu, names_to = "compare", values_to = "p.adj") %>%
    separate(col = compare, c("group1", "group2"), sep = "-", remove = TRUE)
}
#Function to tidy pairwise comparisons to use for adding stats to genus plots----
tidy_pairwise_genus <- function(spread_pairwise){
  spread_pairwise %>%
    pivot_longer(-genus, names_to = "compare", values_to = "p.adj") %>%
    separate(col = compare, c("group1", "group2"), sep = "-", remove = TRUE)
}
#Function to calculate the median shannon values from a dataframe (x) grouped by treatment
get_shannon_median_group <- function(x){
  x %>%
    group_by(group) %>%
    summarize(median=median(shannon)) %>%
    spread(key=group, value=median)
}

#Function to calculate the median sobs (richness) values from a dataframe (x) grouped by treatment
get_sobs_median_group <- function(x){
  x %>%
    group_by(group) %>%
    summarize(median=median(sobs)) %>%
    spread(key=group, value=median)
}

#Function to calculate the median agg_rel_abund values from a dataframe (x) grouped by treatment
get_rel_abund_median_group <- function(x){
  x %>%
    group_by(group) %>%
    summarize(median=median(agg_rel_abund)) %>%
    spread(key=group, value=median)
}

#Function to calculate the median agg_rel_abund values from a dataframe (x) grouped by day
get_rel_abund_median_day <- function(x){
  x %>%
    group_by(day) %>%
    summarize(median=median(agg_rel_abund)) %>%
    spread(key=day, value=median)
}

get_rel_abund_median_clearance_status <- function(x){
  x %>%
    group_by(clearance_status_d10) %>%
    summarize(median=median(agg_rel_abund)) %>%
    spread(key=clearance_status_d10, value=median)
}

#Function to pull significant days (adjusted p value < 0.05) after Kruskal-Wallis statistical analysis
pull_sig_days <- function(dataframe){
  dataframe %>% 
    filter(p.value.adj <= 0.05) %>%
    pull(day)
}  
#Function to pull significant taxa (adjusted p value < 0.05) after statistical analysis
pull_significant_taxa <- function(dataframe, taxonomic_level){
  dataframe %>%
    filter(p.value.adj <= 0.05) %>%
    pull({{ taxonomic_level }}) #Embracing transforms taxonomic_level argument into a column name
}

#Function to tidy pairwise histology comparisons to use for adding stats to plots
tidy_histology_pairwise <- function(spread_pairwise){
  spread_pairwise %>%
    pivot_longer(-Tissue, names_to = "compare", values_to = "p.adj") %>%
    separate(col = compare, c("group1", "group2"), sep = "-", remove = TRUE)
}

#Functions related to 16S rRNA sequencing analysis----
#Function to format distance matrix generated with mothur for use in R.
#Source: Sze et al. mSphere 2019 https://github.com/SchlossLab/Sze_PCRSeqEffects_mSphere_2019/blob/master/code/vegan_analysis.R
read_dist <- function(dist_file_name){
  
  linear_data <- scan(dist_file_name, what="character", sep="\n", quiet=TRUE)
  
  n_samples <- as.numeric(linear_data[1])
  linear_data <- linear_data[-1]
  
  samples <- str_replace(linear_data, "\t.*", "")
  linear_data <- str_replace(linear_data, "[^\t]*\t", "")
  linear_data <- linear_data[-1]
  
  distance_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
  
  for(i in 1:(n_samples-1)){
    row <- as.numeric(unlist(str_split(linear_data[i], "\t")))
    distance_matrix[i+1,1:length(row)] <- row
  }
  
  distance_matrix <- distance_matrix + t(distance_matrix)
  rownames(distance_matrix) <- samples
  
  as.dist(distance_matrix)
}

#Function to find which significant otus/genera/families are shared
intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}

#Shape scheme for differntiating types of tissues----
sample_shape_scheme <- c(19, 17, 15)
sample_shape_groups <- c("cecum", "proximal_colon", "distal_colon")
sample_shape_labels <- c("cecum", "proximal colon", "distal colon")

#Functions to plot data----
#Function to plot cfu data
#Arguments: df = dataframe to plot
#When using the function can add line to specify x axis scale (scale_x_continuous())
plot_cfu_data <- function(df){
  median_summary <- df %>%
    group_by(group, day) %>%
    summarize(median_avg_cfu = median(avg_cfu, na.rm = TRUE))
  #Plot cfu for just the inital 10days
  cfu_plot <- ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = avg_cfu, color= group, fill = group, alpha = day), size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(median_summary, mapping = aes(x = day, y = median_avg_cfu, group = group, color = group), alpha = 0.6, size = 1.5) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    theme_classic()+
    labs(x = "Days Post-Infection", y = "CFU/g Feces") +
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12))+ #Scientific notation labels for y-axis
    geom_hline(yintercept = 100, linetype=2) + #Line that represents our limit of detection when quantifying C. difficile CFU by plating
    geom_text(x = 11, y = 104, color = "black", label = "LOD") + #Label for line that represents our limit of detection when quantifying C. difficile CFU by plating
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme(legend.position = "bottom",
          legend.key= element_rect(colour = "transparent", fill = "transparent"),
          text = element_text(size = 16), # Change font size for entire plot
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days)
}

#Function to plot weight.
#Arguments: df = dataframe you want to plot
#When using the function can add line to specify x axis scale (scale_x_continuous())
plot_weight <- function(df){
  median_summary <- df %>%
    group_by(group, day) %>%
    summarize(median_weight_change = median(weight_change, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = weight_change, color= group, fill = group, alpha = day), size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(median_summary, mapping = aes(x = day, y = median_weight_change, color = group), alpha = 0.6, size = 1.5) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    labs(x = "Days Post-Infection", y = "Weight Change (g)") +
    ylim(-6, 4)+ #Make y-axis for weight_change data uniform across figures
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme_classic()+
    theme(legend.position = "none", #Get rid of legend
          axis.ticks.x = element_blank(),
          text = element_text(size = 16), # Change font size for entire plot
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days)
}

#Simplified function that only plots the median line for each group.
plot_weight_medians <- function(df){
  median_summary <- df %>%
    group_by(group, day) %>%
    summarize(median_weight_change = median(weight_change, na.rm = TRUE))
  ggplot(NULL) +
    geom_line(median_summary, mapping = aes(x = day, y = median_weight_change, color = group), alpha = 0.6, size = 1.5) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    labs(x = "Days Post-Infection", y = "Weight Change (g)") +
    ylim(-6, 4)+ #Make y-axis for weight_change data uniform across figures
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme_classic()+
    theme(legend.position = "none", #Get rid of legend
          axis.ticks.x = element_blank(),
          text = element_text(size = 16), # Change font size for entire plot
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days)
}

#Function to Plot Shannon Diversity Overtime
plot_shannon_overtime <- function(df) {
  median_summary <- df %>%
    group_by(group, day) %>%
    mutate(median_shannon = median(shannon)) %>%
    ggplot(x = day, y = shannon, colour = group)+
    geom_point(mapping = aes(x = day, y = shannon, group = group, color = group, fill = group, alpha = day), size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mapping = aes(x = day, y = median_shannon, group = group, color = group), alpha = 0.6, size = 1) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels) +
    scale_y_continuous(limits = c(0,4.1), expand = c(0, .2))+ #expand argument gets rid of the extra space around the scale
    theme_classic()+
    labs(title=NULL,
         x="Days Post-Infection",
         y="Shannon Diversity Index")+
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme(legend.position = "bottom",
          text = element_text(size = 16), # Change font size for entire plot
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days))
}

#Function to Plot Shannon Diversity Overtime in tissue samples
plot_shannon_overtime_t <- function(df) {
  median_summary <- df %>%
    group_by(group, day) %>%
    mutate(median_shannon = median(shannon)) %>%
    ggplot(x = day, y = shannon, colour = group)+
    geom_point(mapping = aes(x = day, y = shannon, group = group, color = group, fill = group, alpha = day, shape = sample_type), size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mapping = aes(x = day, y = median_shannon, group = group, color = group), alpha = 0.6, size = 1) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels) +
    scale_shape_manual(values=sample_shape_scheme,
                       breaks=sample_shape_groups,
                       labels=sample_shape_labels)+
    scale_y_continuous(limits = c(0,4.1), expand = c(0, .2))+ #expand argument gets rid of the extra space around the scale
    theme_classic()+
    labs(title=NULL,
         x="Days Post-Infection",
         y="Shannon Diversity Index")+
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme(legend.position = "bottom",
          text = element_text(size = 16), # Change font size for entire plot
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days))
}

#Function to Plot Richness (sobs) Overtime
plot_richness_overtime <- function(df) {
  median_summary <- df %>%
    group_by(group, day) %>%
    mutate(median_sobs = median(sobs)) %>%
    ggplot(x = day, y = sobs, group = group, colour = group)+
    geom_point(mapping = aes(x = day, y = sobs, group = group, color = group, fill = group, alpha = day), size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mapping = aes(x = day, y = median_sobs, group = group, color = group), alpha = 0.6, size = 1) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels) +
    scale_y_continuous(limits = c(0,160), expand = c(0, .2))+ #expand argument gets rid of the extra space around the scale
    theme_classic()+
    labs(title=NULL,
         x="Days Post-Infection",
         y="Number of Observed OTUs")+
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme(legend.position = "bottom",
          text = element_text(size = 16), # Change font size for entire plot
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days))
}

#Function to Plot Richness (sobs) Overtime in tissue samples
plot_richness_overtime_t <- function(df) {
  median_summary <- df %>%
    group_by(group, day) %>%
    mutate(median_sobs = median(sobs)) %>%
    ggplot(x = day, y = sobs, group = group, colour = group)+
    geom_point(mapping = aes(x = day, y = sobs, group = group, color = group, fill = group, alpha = day, shape = sample_type), size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mapping = aes(x = day, y = median_sobs, group = group, color = group), alpha = 0.6, size = 1) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels) +
    scale_shape_manual(values=sample_shape_scheme,
                       breaks=sample_shape_groups,
                       labels=sample_shape_labels)+
    scale_y_continuous(limits = c(0,160), expand = c(0, .2))+ #expand argument gets rid of the extra space around the scale
    theme_classic()+
    labs(title=NULL,
         x="Days Post-Infection",
         y="Number of Observed OTUs")+
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme(legend.position = "bottom",
          text = element_text(size = 16), # Change font size for entire plot
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days))
}

#Function to plot PCoA data
plot_pcoa <- function(df){
  ggplot(df, aes(x=axis1, y=axis2, color = group, alpha = day)) +
    geom_point(size=2) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels, 
                        guide = "none")+ #Suppress legend for group colors
    coord_fixed() +
    labs(x="PCoA 1",
         y="PCoA 2",
         alpha= "Day") +
    theme_classic()+
    theme(legend.position = "bottom",
          text = element_text(size = 16))
}

#Function to plot PCoA data from tissue samples, differentiate sample_type with shapes
plot_pcoa_t <- function(df){
  ggplot(df, aes(x=axis1, y=axis2, color = group, alpha = day, shape = sample_type)) +
    geom_point(size=2) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels, 
                        guide = "none")+ #Suppress legend for group colors
    scale_shape_manual(values=sample_shape_scheme,
                       breaks=sample_shape_groups,
                       labels=sample_shape_labels)+
    coord_fixed() +
    labs(x="PCoA 1",
         y="PCoA 2",
         alpha= "Day") +
    theme_classic()+
    theme(legend.position = "bottom",
          text = element_text(size = 16))
}

#Function to format Kruskal-Wallis data frame with adjusted p values to use as a label on plots of cfu, weight, diversity, and other variables over time:
kw_label <- function(dataframe){
  dataframe %>% 
    filter(p.value.adj <= 0.05) %>%
    mutate(p.signif = case_when(
      p.value.adj > 0.05 ~ "NS",
      p.value.adj <= 0.05 ~ "*"
    )) %>%
    pull(p.signif)
}

#Function to plot a list of OTUs across groups of mice at a specific timepoint:
#Arguments: 
#sample_df = subset dataframe of samples to be plotted
#otus = list of otus to plot
#timepoint = day of the experiment to plot
plot_otus_dx <- function(sample_df, otus, timepoint){
  sample_df %>%
    filter(otu %in% otus) %>%
    filter(day == timepoint) %>%
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% # 2,000 is 2 times the subsampling parameter of 1000
    ggplot(aes(x= otu_name, y=agg_rel_abund, color=group))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    geom_hline(yintercept=1/1000, color="gray")+
    stat_summary(fun = 'median',
                 fun.max = function(x) quantile(x, 0.75),
                 fun.min = function(x) quantile(x, 0.25),
                 position = position_dodge(width = 1)) +
    labs(title=NULL,
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "none",
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Function to plot a list of genera across groups of mice at a specific timepoint:
#Arguments: 
#sample_df = subset dataframe of samples to be plotted
#genus = list of genera to plot
#timepoint = day of the experiment to plot
plot_genus_dx <- function(sample_df, genera, timepoint){
  sample_df %>%
    filter(genus %in% genera) %>%
    filter(day == timepoint) %>%
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% # 2,000 is 2 times the subsampling parameter of 1000
    ggplot(aes(x= genus, y=agg_rel_abund, color=group))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    geom_hline(yintercept=1/1000, color="gray")+ #Represents limit of detection
    stat_summary(fun = 'median',
                 fun.max = function(x) quantile(x, 0.75),
                 fun.min = function(x) quantile(x, 0.25),
                 position = position_dodge(width = 1)) +
    labs(title=NULL,
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "none",
          axis.text.y = element_text(face = "italic"), #Have only the genus names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Function to create a heatmap plot the relative abundances of a list of OTUs over time, faceted by group----
#Arguments: 
#sample_df = subset dataframe of samples to be plotted
#otus = list of otus to plot
#timepoints = days of the experiment to plot
hm_plot_otus <- function(sample_df, otus, timepoints){
  sample_df %>%
    mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
    mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                             "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the groups  
    filter(otu %in% otus) %>%
    filter(day %in% timepoints) %>% 
    group_by(group, otu_name, day) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_tile(aes(x = day, y=otu_name, fill=median))+
    labs(title=NULL,
         x=NULL,
         y=NULL)+
    facet_wrap(~group, labeller = labeller(group = facet_labels)) + #Make sure you specify facet_labels before running function
    #    scale_fill_gradient2(low="white", mid=color_scheme, high = 'black',
    #                         limits = c(1/10000, 1), na.value = NA, midpoint = .3,
    #                         breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) + 
    scale_fill_distiller(trans = "log10",palette = "YlGnBu", direction = 1, name = "Relative \nAbundance",
                         limits = c(1/10000, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    #    scale_y_discrete(limits=rev(levels(as.factor(sample_df$otu_name))))+#List OTU names alphabetically
    theme(plot.title=element_text(hjust=0.5),
          strip.background = element_blank(), #get rid of box around facet_wrap labels
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Function to create a heatmap plot the relative abundances of a list of Genera over time, faceted by group----
#Arguments: 
#sample_df = subset dataframe of samples to be plotted
#genera = list of genera to plot
#timepoints = days of the experiment to plot
hm_plot_genus <- function(sample_df, genera, timepoints){
  sample_df %>%
    mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
    mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                             "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the groups  
    filter(genus %in% genera) %>%
    filter(day %in% timepoints) %>% 
    group_by(group, genus, day) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_tile(aes(x = day, y=genus, fill=median))+
    labs(title=NULL,
         x=NULL,
         y=NULL)+
    facet_wrap(~group, labeller = labeller(group = facet_labels)) + #Make sure you specify facet_labels before running function
    #    scale_fill_gradient2(low="white", mid=color_scheme, high = 'black',
    #                         limits = c(1/10000, 1), na.value = NA, midpoint = .3,
    #                         breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) + 
    scale_fill_distiller(trans = "log10",palette = "YlGnBu", direction = 1, name = "Relative \nAbundance",
                         limits = c(1/10000, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    #    scale_y_discrete(limits=rev(levels(as.factor(sample_df$genus))))+#List genera names alphabetically
    theme(plot.title=element_text(hjust=0.5),
          strip.background = element_blank(), #get rid of box around facet_wrap labels
          axis.text.y = element_markdown(face = "italic"), #Have only the genus names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}
#Function to create a heatmap plot the relative abundances of a list of Genera over time, faceted by genus----
#Arguments: 
#sample_df = subset dataframe of samples to be plotted
#exp_groups = character vector of group abbr. 
#Genera = list of genera to plot
#timepoints = days of the experiment to plot
#exp_group_labels = character vector of full group name 
hm_plot_genus_facet <- function(sample_df, exp_groups, genera_list, timepoints, exp_group_labels){
  sample_df %>%
    mutate(group = fct_relevel(group, exp_groups)) %>% #Specify the order of the groups
    mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
    mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                             "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the days  
    filter(genus %in% genera_list) %>%
    filter(day %in% timepoints) %>% 
    group_by(group, genus, day) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_tile(aes(x = day, y=group, fill=median))+
    labs(title=NULL,
         x=NULL,
         y=NULL)+
    facet_wrap(~genus, nrow = 2, labeller = label_wrap_gen(width = 10)) + 
    scale_fill_distiller(trans = "log10",palette = "YlGnBu", direction = 1, name = "Relative \nAbundance",
                         limits = c(1/10000, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    scale_y_discrete(label = exp_group_labels)+ #Descriptive group names that match the rest of the plots
    theme_classic()+
    theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
          strip.text = element_text(face = "italic"),
          plot.title = element_markdown(hjust = 0.5), #Have only the genera names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Function to create a heatmap plot the relative abundances of a list of OTUs over time, faceted by group----
#Arguments: 
#sample_df = subset dataframe of samples to be plotted
#otus = list of otus to plot
#timepoints = days of the experiment to plot
hm_plot_tissues <- function(sample_df, otus, timepoints){
  sample_df %>%
    mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
    mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                             "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the groups  
    filter(otu %in% otus) %>%
    filter(day %in% timepoints) %>% 
    group_by(group, otu_name, day, sample_type) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_tile(aes(x = day, y=otu_name, fill=median))+
    labs(title=NULL,
         x=NULL,
         y=NULL)+
    facet_wrap(~sample_type, labeller = labeller(sample_type = facet_labels)) + #Make sure you specify facet_labels before running function
    scale_fill_distiller(trans = "log10",palette = "YlGnBu", direction = 1, name = "Relative \nAbundance",
                         limits = c(1/10000, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          strip.background = element_blank(), #get rid of box around facet_wrap labels
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Function to create a heatmap plot the relative abundances of a list of genera over time, faceted by group----
#Arguments: 
#sample_df = subset dataframe of samples to be plotted
#genera = list of genera to plot
#timepoints = days of the experiment to plot
hm_plot_tissues_genera <- function(sample_df, genera, timepoints){
  sample_df %>%
    mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
    mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                             "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the groups  
    filter(genus %in% genera) %>%
    filter(day %in% timepoints) %>% 
    group_by(group, genus, day, sample_type) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_tile(aes(x = day, y=genus, fill=median))+
    labs(title=NULL,
         x=NULL,
         y=NULL)+
    facet_wrap(~sample_type, labeller = labeller(sample_type = facet_labels)) + #Make sure you specify facet_labels before running function
    scale_fill_distiller(trans = "log10",palette = "YlGnBu", direction = 1, name = "Relative \nAbundance",
                         limits = c(1/10000, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          strip.background = element_blank(), #get rid of box around facet_wrap labels
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Function to create a heatmap demonstrating the relative abundance of 1 OTU over time across all groups----
#sample_df = subset dataframe of samples to be plotted
#specify_otu = 1 OTU to plot (name should be in quotes)
#timepoints = days of the experiment to plot
hm_1_otu <- function(sample_df, specify_otu, timepoints){
  otu_format_name <-sample_df %>% select(otu, otu_name) %>% distinct(otu, otu_name) %>% 
    filter(otu == specify_otu) %>% pull(otu_name) #Get correctly formatted otu name to use with element_markdown
  sample_df %>%
    mutate(group = fct_relevel(group, "CN", "C", "1RM1", "M1", "FRM", "RM", "CWM", "WMR", "WMC", "WMN", "WM")) %>% #Specify the order of the groups
    mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
    mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                             "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the groups  
    filter(otu == specify_otu) %>%
    filter(day %in% timepoints) %>% 
    group_by(group, day) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_tile(aes(x = day, y=group, fill=median))+
  labs(title=otu_format_name,
       x=NULL,
       y=NULL)+
  scale_fill_distiller(trans = "log10",palette = "YlGnBu", direction = 1, name = "Relative \nAbundance",
                       limits = c(1/10000, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  scale_y_discrete(label = c("Clind. without infection", "Clind.", "1-day PEG 3350 + 1-day recovery", "1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350",
                             "Clind. + 1-day PEG 3350", "5-day PEG 3350 + 10-day recovery", "5-day PEG 3350 + Clind.", 
                             "5-day PEG 3350 without infection", "5-day PEG 3350"))+ #Descriptive group names that match the rest of the plots
  theme_classic()+
  theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
        plot.title = element_markdown(hjust = 0.5), #Have only the OTU names show up as italics
        text = element_text(size = 16)) # Change font size for entire plot
}

#Function to create a heatmap demonstrating the relative abundance of 1 genus over time across all groups----
#sample_df = subset dataframe of samples to be plotted
#specify_genus = 1 genus to plot (name should be in quotes)
#timepoints = days of the experiment to plot
hm_1_genus <- function(sample_df, specify_genus, timepoints){
  genus_format_name <-sample_df %>% select(genus) %>% distinct(genus) %>% 
    filter(genus == specify_genus) %>% pull(genus) #Get genus name
  sample_df %>%
    mutate(group = fct_relevel(group, "CN", "C", "1RM1", "M1", "FRM", "RM", "CWM", "WMR", "WMC", "WMN", "WM")) %>% #Specify the order of the groups
    mutate(day = factor(day, levels = unique(as.factor(day)))) %>% #Transform day variable into factor variable
    mutate(day = fct_relevel(day, "-15", "-11", "-10", "-5", "-4", "-2", "-1", "0", "1", "2", "3", "4",
                             "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% #Specify the order of the groups  
    filter(genus == specify_genus) %>%
    filter(day %in% timepoints) %>% 
    group_by(group, day) %>% 
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>%  #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_tile(aes(x = day, y=group, fill=median))+
    labs(title=genus_format_name,
         x=NULL,
         y=NULL)+
    scale_fill_distiller(trans = "log10",palette = "YlGnBu", direction = 1, name = "Relative \nAbundance",
                         limits = c(1/10000, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    scale_y_discrete(label = c("Clind. without infection", "Clind.", "1-day PEG 3350 + 1-day recovery", "1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350",
                               "Clind. + 1-day PEG 3350", "5-day PEG 3350 + 10-day recovery", "5-day PEG 3350 + Clind.", 
                               "5-day PEG 3350 without infection", "5-day PEG 3350"))+ #Descriptive group names that match the rest of the plots
    theme_classic()+
    theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
          plot.title = element_text(face = "italic", hjust = 0.5), #Italicize genus name
          text = element_text(size = 16)) # Change font size for entire plot
}

#Function to plot a genus over time
#genus_plot = genus to plot in quotes. Ex: "Bacteroides"
#sample_df = subset dataframe of just stool or tissue samples
genus_over_time <- function(genus_plot, sample_df){
  specify_genus_name <- sample_df %>% 
    filter(genus == genus_plot) %>% 
    pull(genus)
  genus_median <- sample_df %>% 
    filter(genus == genus_plot) %>% 
    group_by(group, day) %>% 
    summarize(median=(median(agg_rel_abund + 1/2000))) %>% 
    ungroup
  genus_mice <- sample_df %>% 
    filter(genus == genus_plot) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>%
    select(day, agg_rel_abund, genus, group)
  genus_time <- ggplot(NULL)+
    geom_point(genus_mice, mapping = aes(x=day, y=agg_rel_abund, color=group), size  = 1.5, position = position_dodge(width = 0.6))+
    geom_line(genus_median, mapping = aes(x=day, y=median, group = group, color=group), size = 1, show.legend = FALSE)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    geom_hline(yintercept=1/1000, color="gray")+
    labs(title=specify_genus_name,
         x="Day",
         y="Relative abundance (%)") +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_text(hjust = 0.5, face = "italic"),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"),  # Add gray lines to clearly separate symbols by days)
          text = element_text(size = 18)) # Change font size for entire plot
}

#Function to plot an otu_over_time
#otu_plot = otu to plot in quotes. Ex: "Peptostreptococcaceae (OTU 12)"
#sample_df = subset dataframe of just stool or tissue samples
otu_over_time <- function(otu_plot, sample_df){
  specify_otu_name <- sample_df %>% 
    filter(otu == otu_plot) %>% 
    pull(otu_name)
  otu_median <- sample_df %>% 
    filter(otu == otu_plot) %>% 
    group_by(group, day) %>% 
    summarize(median=(median(agg_rel_abund + 1/2000))) %>% 
    ungroup
  otu_mice <- sample_df %>% 
    filter(otu == otu_plot) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>%
    select(day, agg_rel_abund, otu, group)
  otu_time <- ggplot(NULL)+
    geom_point(otu_mice, mapping = aes(x=day, y=agg_rel_abund, color=group), size  = 1.5, position = position_dodge(width = 0.6))+
    geom_line(otu_median, mapping = aes(x=day, y=median, group = group, color=group), size = 1, show.legend = FALSE)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    geom_hline(yintercept=1/1000, color="gray")+
    labs(title=specify_otu_name,
         x="Day",
         y="Relative abundance (%)") +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_markdown(hjust = 0.5),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"),  # Add gray lines to clearly separate symbols by days)
          text = element_text(size = 18)) # Change font size for entire plot
}

#Examine a specific OTU at a single timepoint across different sample types in one group of mice
#otu_plot = name of otu to plot
#sample_df = dataframe of samples to test
#timepoint = experiment day to examine "30"
#group_name = group to examine, "WMR"
otu_gi_distrib <- function(otu_plot, sample_df, timepoint, group_name){
  specify_otu_name <- sample_df %>% 
    filter(otu == otu_plot) %>% 
    pull(otu_name)
  otu_median <- sample_df %>% 
    mutate(mouse_id = factor(mouse_id, levels = unique(as.factor(mouse_id)))) %>% 
    filter(otu == otu_plot & day == timepoint, group == group_name) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>%
    group_by(sample_type) %>% 
    mutate(median=(median(agg_rel_abund + 1/2000))) %>% 
    ungroup() %>% 
    ggplot(aes(x=sample_type, y = agg_rel_abund, color = mouse_id))+    
    geom_jitter(shape = 19, size=2, alpha = 0.5) +    
    geom_errorbar(aes(ymax= median, ymin= median), color = "gray50", size = 1)+#Add line to show median of each sample type
    geom_hline(yintercept=1/1000, color="gray")+
    labs(title=specify_otu_name,
         color = "Mouse ID",
         x=NULL,
         y="Relative abundance (%)") +
    scale_x_discrete(label = c("Cecum", "Proximal colon", "Distal colon"))+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_markdown(hjust = 0.5),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"),  # Add gray lines to clearly separate symbols by days)
          text = element_text(size = 18)) # Change font size for entire plot
}

#Function to create faceted line plots of relative abundances of genera of interest over time
#sample_df = subset dataframe of samples to be plotted
#genera = list of genera to plot
#timepoints = days of the experiment to plot
#specify_linetype = "solid" for C. diff challenged. "dashed" for mice that were mock challenged
line_plot_genus <- function(sample_df, genera, timepoints, specify_linetype){
  sample_df %>% 
    filter(genus %in% genera) %>% #Select only genera of interest
    mutate(genus = fct_relevel(genus, genera)) %>% #Reorder genera to match order of genera of interest
    mutate(group = fct_relevel(group, rev(color_groups))) %>% #Specify the order of the groups
    filter(day %in% timepoints) %>% #Select only timepoints of interest
    group_by(group, genus, day) %>% 
    mutate(day = as.integer(day)) %>%
    summarize(median=median(agg_rel_abund + 1/2000),`.groups` = "drop") %>% #Add small value (1/2Xsubssampling parameter) so that there are no infinite values with log transformation
    ggplot()+
    geom_line(aes(x = day, y=median, color=group), linetype = specify_linetype)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    scale_x_continuous(limits = c(-1.5,11), breaks = c(-1:10), labels = c(-1:10))+
    scale_y_continuous(trans = "log10", limits = c(1/10900, 1), breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    geom_hline(yintercept=1/1000, color="gray")+ #Represents limit of detection
    labs(title=NULL,
         x="Days Post-Infection",
         y="Relative abundance (%)")+
    facet_wrap(~genus, nrow = 2, labeller = label_wrap_gen(width = 10))+
    theme_classic()+
    theme(strip.background = element_blank(), #get rid of box around facet_wrap labels
          strip.text = element_text(face = "italic"),
          plot.title = element_markdown(hjust = 0.5), #Have only the genera names show up as italics
          text = element_text(size = 16),
          legend.position = "None")
}
