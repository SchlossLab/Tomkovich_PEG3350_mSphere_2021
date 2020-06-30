#Load required libraries
library(tidyverse)
library(readxl)
library(writexl)
library(cowplot)
library(broom)

#Read in metadata
metadata <- read_excel("data/process/metadata.xlsx", col_types = c("numeric", "text", "numeric", "text", "numeric", "text", "text", "text", "text", "text", "numeric", "numeric")) %>% #specify column types
  rename(mouse_id = `Mouse ID`,
         group = Group,
         cage = `Cage #`,
         ear_mark = `Ear Mark`,
         day = Day,
         weight = Weight,
         avg_cfu = avgCFU) %>% #get rid of spaces in column names and make all variables lowercase
  unite(col = m_id_unique, c(mouse_id, exp_num), sep = "_", remove = FALSE) #add column to differentiate individual mice by merging mouse ID # and experiment number

#Check for duplicated unique_labels
duplicated <- metadata %>% 
  filter(duplicated(unique_label)) #0 duplicates

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

#Function to tidy pairwise comparisons to use for adding stats to plots----
tidy_pairwise <- function(spread_pairwise){
  spread_pairwise %>%
    pivot_longer(-day, names_to = "compare", values_to = "p.adj") %>%
    separate(col = compare, c("group1", "group2"), sep = "-", remove = TRUE)
}
