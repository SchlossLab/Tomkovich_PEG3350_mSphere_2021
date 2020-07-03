#Load required libraries
library(tidyverse)
library(readxl)
library(writexl)
library(cowplot)
library(broom)
library(ggpubr)

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

#Function to tidy pairwise histology comparisons to use for adding stats to plots
tidy_histology_pairwise <- function(spread_pairwise){
  spread_pairwise %>%
    pivot_longer(-Tissue, names_to = "compare", values_to = "p.adj") %>%
    separate(col = compare, c("group1", "group2"), sep = "-", remove = TRUE)
}

#Functions used to plot data----
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
    theme_classic()
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
    theme_classic()
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
    theme_classic()
}
