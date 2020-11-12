#Load required libraries
library(tidyverse)
library(readxl)
library(writexl)
library(cowplot)
library(broom)
library(ggpubr)
library(gganimate)
library(glue)
library(ggtext)
library(vegan)

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
  anti_join(duplicates_to_drop)

#Join metadata to 16S seq prep metadata----
metadata <- seq_prep_metadata %>% 
  select(unique_label, ext_plate, miseq_run) %>% #Just select unique_label and the plate # and miseq run variables to test with adonis/PERMANOVA
  full_join(metadata, seq_prep_metadata, by = "unique_label") %>% 
  mutate(sample_type = case_when(str_detect(unique_label, "water") ~ "water",
                                 str_detect(unique_label, "FMT") ~ "FMT",
                                 TRUE ~ sample_type),
         group = case_when(str_detect(unique_label, "water") ~ "water",
                           str_detect(unique_label, "FMT") ~ "FMT",
                           TRUE ~ group)) %>% 
  mutate(group = factor(group, levels = unique(as.factor(group)))) #Transform group variable into factor variable

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
    geom_point(df, mapping = aes(x = day, y = avg_cfu, color= group, fill = group), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(median_summary, mapping = aes(x = day, y = median_avg_cfu, color = group), alpha = 0.6, size = 1.5) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    labs(x = "Days Post-Infection", y = "CFU/g Feces") +
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12))+ #Scientific notation labels for y-axis
    geom_hline(yintercept = 100, linetype=2) + #Line that represents our limit of detection when quantifying C. difficile CFU by plating
    geom_text(x = 11, y = 104, color = "black", label = "LOD") + #Label for line that represents our limit of detection when quantifying C. difficile CFU by plating
    theme(text = element_text(size = 16))+  # Change font size for entire plot
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme_classic()+
    theme(legend.position = "bottom",
          legend.key= element_rect(colour = "transparent", fill = "transparent"),
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
    geom_point(df, mapping = aes(x = day, y = weight_change, color= group, fill = group), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(median_summary, mapping = aes(x = day, y = median_weight_change, color = group), alpha = 0.6, size = 1.5) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    labs(x = "Days Post-Infection", y = "Weight Change (g)") +
    ylim(-6, 4)+ #Make y-axis for weight_change data uniform across figures
    theme(text = element_text(size = 16))+  # Change font size for entire plot
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme_classic()+
    theme(legend.position = "none", #Get rid of legend
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
    theme(text = element_text(size = 16))+  # Change font size for entire plot
    annotate("text", y = y_position, x = x_annotation, label = label, size =7)+ #Add statistical annotations
    theme_classic()+
    theme(legend.position = "none", #Get rid of legend
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days)
}

#Function to Plot Shannon Diversity Overtime
plot_shannon_overtime <- function(df) {
  median_summary <- df %>%
    group_by(group, day) %>%
    mutate(median_shannon = median(shannon)) %>%
    ggplot(x = day, y = shannon, colour = group)+
    geom_point(mapping = aes(x = day, y = shannon, color = group, fill = group), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mapping = aes(x = day, y = median_shannon, color = group), alpha = 0.6, size = 1) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels) +
    theme_classic()+
    theme(legend.position = c(.9,.25),
          text = element_text(size = 14), # Change font size for entire plot
          axis.ticks.x = element_blank())
}

#Function to plot PCoA data
plot_pcoa <- function(df){
  ggplot(df, aes(x=axis1, y=axis2, color = group, alpha = day)) +
    geom_point(size=2) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    #    scale_alpha_continuous(range = c(.3, 1),
    #                           breaks= c(2, 4, 6, 8, 10),
    #                           labels=c(2, 4, 6, 8, 10))+
    coord_fixed() +
    #    xlim(-0.4, 0.65)+
    #    ylim(-0.45, 0.6)+
    labs(x="PCoA 1",
         y="PCoA 2",
         alpha= "Day") +
    theme_classic()
}




