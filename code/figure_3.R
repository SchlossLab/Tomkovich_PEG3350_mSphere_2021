source("code/functions.R") #Reads in metadata

#Define color scheme----
color_scheme <- c("#225ea8", "#238b45", "#88419d") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c( "1RM1", "C", "M1")
color_labels <- c( "1-day PEG 3350 + 1-day recovery", "Clind.", "1-day PEG 3350")

#Narrow metadata to relevant groups and experiments (C, 1RM1, M1)----
fig3_metadata <- metadata %>% 
  filter(group == "C" & exp_num %in% c("M6")| #Only use C mice from this experiments. Allocated groups to figures based on paper outline.
         group == "1RM1" & exp_num %in% c("M6R")| #Had to differentiate experiment 6 from 6R in the metadata to create m_id_unique that wouldn't overlap for the M1 & 1RM1 mice that are both labeled with mouse_ids that are #s1-6
         group == "M1" & exp_num %in% c("M6"))%>% 
  mutate(group=factor(group, levels=c("C", "1RM1", "M1")))  # Make sure group is treated as a factor


# of mice represented in the figure
fig3_mice <- length(unique(fig3_metadata$m_id_unique)) 
# 18 mice total for figure 3

#Note baseline weight for each group of mice (based on the earliest timepoint recorded for each experiment)----
baseline <- fig3_metadata %>% #Baseline weight was taken at day -5 for groups C, WM, and WMC
  filter(group == "C" & day == -15| #6 mice in C group
         group == "1RM1" & day == -2| #6 mice in 1RM1 group
         group == "M1" & day == -11) %>%  #6 mice in M1 group
  mutate(baseline_weight = weight) %>% #This column represents the initial weight that was recorded for each mouse
  select(m_id_unique, baseline_weight) #Will use m_id_unique to join baseline_weights to fig3_metadata

#Make a new column that represents weight_change from baseline_weight----
fig3_weightdata <- inner_join(fig3_metadata, baseline, by = "m_id_unique") %>% #Join baseline weight to fig3_metadata
  group_by(m_id_unique, day) %>% #Group by each unique mouse and experiment day
  mutate(weight_change = weight-baseline_weight) %>% #Make a new column that represents the change in weight from baseline (all weights recorded in grams)
  ungroup() %>% 
  filter(!is.na(weight)) #drop rows with NA values for fig3_weightdata. 378 samples including NAs, 306 samples after excluding NAs

#C. diff CFU plot----
#Narrow fig3_metadata to just timepoints relevant to C. difficile CFU tracking (Anything on or after day 0)
fig3_cfudata <- fig3_metadata %>% 
  filter(day > -1)
fig3_cfu_na <- sum(is.na(fig3_cfudata$avg_cfu)) #53 samples with NA values. 5 samples, where we weren't able to get a stool sample. Rest of NAs are from timepoints after D15 which is the day we stopped tracking C. diff CFU
#Drop rows with NA values for fig3_cfu:
fig3_cfudata <- fig3_cfudata %>% 
  filter(!is.na(avg_cfu)) #181 samples total

plot_cfu <- function(df){
  mean_summary <- df %>% 
    group_by(group, day) %>% 
    summarize(mean_avg_cfu = mean(avg_cfu, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = avg_cfu, color= group, fill = group), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mean_summary, mapping = aes(x = day, y = mean_avg_cfu, color = group), alpha = 0.6, size = 1) +
    scale_colour_manual(name=NULL,
                   values=color_scheme,
                  breaks=color_groups,
                 labels=color_labels)+
    labs(x = "Days Post-Infection", y = "CFU/g Feces") +
    scale_x_continuous(breaks = c(0, 5, 10),
                       limits = c(-1, 11)) +
    geom_hline(yintercept = 100, linetype=2) + #Line that represents our limit of detection when quantifying C. difficile CFU by plating
    geom_text(x = 11, y = 104, color = "black", label = "LOD") + #Label for line that represents our limit of detection when quantifying C. difficile CFU by plating
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9))+ #Scientific notation labels for y-axis
    theme_classic()
}
fig3_cfu <- plot_cfu(fig3_cfudata)
save_plot(filename = "results/figures/fig3_cfu.png", fig3_cfu, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)


#Weight change plot----
#Function to plot weight. Argument = dataframe you want to plot
plot_weight <- function(df){
  mean_summary <- df %>% 
    group_by(group, day) %>% 
    summarize(mean_weight_change = mean(weight_change, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = weight_change, color= group, fill = group), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mean_summary, mapping = aes(x = day, y = mean_weight_change, color = group), alpha = 0.6, size = 1) +
    scale_colour_manual(name=NULL,
                   values=color_scheme,
                  breaks=color_groups,
                 labels=color_labels)+
    labs(x = "Days Post-Infection", y = "Weight Change (g)") +
    ylim(-6, 4)+ #Make y-axis for weight_change data uniform across figures
    scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10),
                       limits = c(-16, 11)) +
    theme_classic()
}

#Simplified function to plot weight that only plots the mean of each group and no points for individual mice. Argument = dataframe you want to plot.
plot_weight_simple <- function(df){
  mean_summary <- df %>% 
    group_by(group, day) %>% 
    summarize(mean_weight_change = mean(weight_change, na.rm = TRUE))
  ggplot(NULL) +
    geom_line(mean_summary, mapping = aes(x = day, y = mean_weight_change, color = group), alpha = 0.6, size = 1) +
    scale_colour_manual(name=NULL,
                   values=color_scheme,
                  breaks=color_groups,
                 labels=color_labels)+
    labs(x = "Days Post-Infection", y = "Weight Change (g)") +
    ylim(-6, 4)+ #Make y-axis for weight_change data uniform across figures
    scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10),
                       limits = c(-16, 11)) +
    theme_classic()
}

fig3_weight <- plot_weight(fig3_weightdata) 
save_plot(filename = "results/figures/fig3_weight.png", fig3_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
fig3v2_weight <- plot_weight_simple(fig3_weightdata)
save_plot(filename = "results/figures/fig3v2_weight.png", fig3v2_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
kruskal_wallis_cfu <- fig3_cfudata %>% 
  filter(day %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) %>%  #only test days that we have CFU data for 
  group_by(day) %>% 
  do(tidy(kruskal.test(avg_cfu~factor(group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) 
#Timepoints where C. diff CFU is significantly different across the sources of mice
sig_C.diff_CFU_timepoints <- kruskal_wallis_cfu %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day)
