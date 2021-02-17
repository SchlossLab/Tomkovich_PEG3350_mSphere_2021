source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Function to read distance data frame
read_dist_df <- function(dist_file_name){
  linear_data <- scan(dist_file_name, what="character", quiet=TRUE)[-1]
  
  samples <- str_subset(linear_data, "D|F") #Pull out sample ids with D (denotes experiment day) F (is part of the label for the FMT samples)
  n_samples <- length(samples)
  distance_strings <- str_subset(linear_data, "\\.")
  
  distance_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
  colnames(distance_matrix) <- samples
  as_tibble(cbind(rows=samples, distance_matrix)) %>%
    gather(columns, distances, -rows) %>%
    filter(rows < columns) %>%
    arrange(columns, rows) %>%
    mutate(distances = as.numeric(distance_strings))
}

#Function that reads in .dist file and reformats the columns (split day into sep. column), 
#then filters to just the rows that have distances for the same mouse at different time points
#Arguments
# file_path = file path to .dist file
format_dist_df <- function(file_path){
  read_dist_df(file_path) %>% 
    #Remove rows that have FMT samples since we will not include these in the analysis
    filter(!str_detect(rows, "FMT")) %>% 
    filter(!str_detect(columns, "FMT")) %>% 
    #Relabel rows and columns by separating day (D#) label into a separate colony
    separate(rows, into=c("row_m_id", "row_day"), sep = "D",remove = FALSE) %>%  #Separate unique label into unique mouse id & day numbers
    separate(columns, into=c("col_m_id", "col_day"), sep = "D",remove = FALSE) %>%  #Separate unique label into unique mouse id & day numbers
    mutate(row_day = gsub("n", "-", row_day),
           col_day = gsub("n", "-", col_day)) %>% #change all n to - so that the unique_label matches the unique_label in motility_samples
    #Select distances that correspond to samples from the same mouse
    filter(row_m_id == col_m_id) %>% 
    #Join distances to relevant metadata for each mouse
    left_join(select(metadata, unique_label, unique_mouse_id, group, unique_cage_no, exp_num),
              by = c("rows" = "unique_label"))
}

#Function to plot Bray-Curtis Distance relative to baseline for each mouse within each group
#input_df = dataframe of the samples that you want to plot
#timepoint = timepoint of interest relative to baseline (number must be in quotes)
plot_bc_dist <- function(input_df, timepoint){
  input_df %>% 
    filter(row_day %in% c("baseline", timepoint),
           col_day %in% c("baseline", timepoint)) %>% 
    group_by(group) %>% 
    mutate(median = median(distances)) %>% 
    ungroup() %>% 
    ggplot(aes(x = group, y = distances, color = group)) +
    geom_errorbar(aes(ymax=median, ymin=median), color = "gray50", size = 1, show.legend = FALSE)+ #Add lines to indicate median 
    geom_jitter(size=2) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+ 
    labs(title = NULL, x = NULL, y = "Bray-Curtis Distance relative \n to baseline")+
    theme(plot.title = element_text(hjust = 0.5)) +#Center plot title
    theme_classic()+
    lims(y = c(0, 1))+
    theme(legend.position = "none",
          text = element_text(size = 19),# Change font size for entire plot
          axis.title.y = element_text(size = 17),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) #Remove x axis ticks
}

#Format 5-days PEG stool samples distance matrix----
five_d_stools <- format_dist_df("data/process/5_day_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist") %>% 
  #Rename day values to baseline based on the group (WMR group baseline is day = -15, but for the rest of the groups day = -5)
  #Same baseline specified in code/5_days_PEG_cfu_weight.R
  mutate(row_day = case_when(group == "C" & row_day == "-5" ~ "baseline",
                             group == "WM" & row_day == "-5" ~ "baseline",
                             group == "WMC" & row_day == "-5" ~ "baseline",
                             group == "WMR" & row_day == "-15" ~ "baseline",
                             TRUE ~ row_day),
         col_day = case_when(group == "C" & col_day == "-5" ~ "baseline",
                             group == "WM" & col_day == "-5" ~ "baseline",
                             group == "WMC" & col_day == "-5" ~ "baseline",
                             group == "WMR" & col_day == "-15" ~ "baseline",
                             TRUE ~ col_day),
         )

#Define color scheme to match 5_days_PEG plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMR")
color_labels <- c( "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 + 10-day recovery")

#Plots of different timepoints of interest for the 5-days PEG subset----
dn1_baseline <- plot_bc_dist(five_d_stools, "-1")
d0_baseline <- plot_bc_dist(five_d_stools, "0")
d1_baseline <- plot_bc_dist(five_d_stools, "1")
d4_baseline <- plot_bc_dist(five_d_stools, "4")
d6_baseline <- plot_bc_dist(five_d_stools, "6")
d10_baseline <- plot_bc_dist(five_d_stools, "10")
save_plot(filename = "results/figures/5_days_PEG_bc_d10_baseline.png", d10_baseline, base_height = 4, base_width = 4)
d30_baseline <- plot_bc_dist(five_d_stools, "30")
save_plot(filename = "results/figures/5_days_PEG_bc_d30_baseline.png", d30_baseline, base_height = 4, base_width = 4)

#Format 1-day PEG stool samples distance matrix----
one_d_stools <- format_dist_df("data/process/1_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist") %>% 
  #Rename day values to baseline based on the group 
  #Same baseline specified in code/1_day_PEG_cfu_weight.R
  mutate(row_day = case_when(group == "C" & row_day == "-15" ~ "baseline",
                             group == "M1" & row_day == "-11" ~ "baseline",
                             group == "1RM1" & row_day == "-2" ~ "baseline",
                             TRUE ~ row_day),
         col_day = case_when(group == "C" & col_day == "-15" ~ "baseline",
                             group == "M1" & col_day == "-11" ~ "baseline",
                             group == "1RM1" & col_day == "-2" ~ "baseline",
                             TRUE ~ col_day),
  )

#Define color scheme for 1-day PEG subset figures----
color_scheme <- c("#238b45", "#88419d", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "M1", "1RM1")
color_labels <- c("Clind.", "1-day PEG 3350", "1-day PEG 3350 + 1-day recovery")

#Plots of different timepoints of interest for the 1-day PEG subset----
dn1_baseline <- plot_bc_dist(one_d_stools, "-1")
d0_baseline <- plot_bc_dist(one_d_stools, "0")
d1_baseline <- plot_bc_dist(one_d_stools, "1")
d2_baseline <- plot_bc_dist(one_d_stools, "2")
d5_baseline <- plot_bc_dist(one_d_stools, "5")
d7_baseline <- plot_bc_dist(one_d_stools, "7")
save_plot(filename = "results/figures/1_day_PEG_bc_d7_baseline.png", d7_baseline, base_height = 4, base_width = 4)

#Format post-CDI PEG stool samples distance matrix----
post_CDI_stools <- format_dist_df("data/process/post_CDI_PEG/stools/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist") %>% 
  #Rename day values to baseline based on the group (-1 for all groups in this subset is the day before treatment)
  mutate(row_day = case_when(group == "C" & row_day == "-1" ~ "baseline",
                             group == "CWM" & row_day == "-1" ~ "baseline",
                             group == "RM" & row_day == "-1" ~ "baseline",
                             group == "FRM" & row_day == "-1" ~ "baseline",
                             TRUE ~ row_day),
         col_day = case_when(group == "C" & col_day == "-1" ~ "baseline",
                             group == "CWM" & col_day == "-1" ~ "baseline",
                             group == "RM" & col_day == "-1" ~ "baseline",
                             group == "FRM" & col_day == "-1" ~ "baseline",
                             TRUE ~ col_day),
  )

#Define color scheme to match post_CDI_PEG plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "CWM", "FRM", "RM")
color_labels <- c( "Clind.", "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350")

#Plots of different timepoints of interest for the post-CDI PEG subset----
d0_baseline <- plot_bc_dist(post_CDI_stools, "0")
d1_baseline <- plot_bc_dist(post_CDI_stools, "1")
d2_baseline <- plot_bc_dist(post_CDI_stools, "2")
d3_baseline <- plot_bc_dist(post_CDI_stools, "3")
d5_baseline <- plot_bc_dist(post_CDI_stools, "5")
d6_baseline <- plot_bc_dist(post_CDI_stools, "6")
d7_baseline <- plot_bc_dist(post_CDI_stools, "7")
d8_baseline <- plot_bc_dist(post_CDI_stools, "8")
d9_baseline <- plot_bc_dist(post_CDI_stools, "9")
d10_baseline <- plot_bc_dist(post_CDI_stools, "10")
save_plot(filename = "results/figures/post_CDI_PEG_bc_d10_baseline.png", d10_baseline, base_height = 4, base_width = 4)
d15_baseline <- plot_bc_dist(post_CDI_stools, "15")
save_plot(filename = "results/figures/post_CDI_PEG_bc_d15_baseline.png", d15_baseline, base_height = 4, base_width = 4)
d30_baseline <- plot_bc_dist(post_CDI_stools, "30")
save_plot(filename = "results/figures/post_CDI_PEG_bc_d30_baseline.png", d30_baseline, base_height = 4, base_width = 4)


