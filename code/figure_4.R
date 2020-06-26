source("code/functions.R") #Reads in metadata

#Define color scheme----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "CWM", "FRM", "RM")
color_labels <- c( "Clind.", "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350")

#Narrow metadata to relevant groups and experiments (C, CWM, RM, FRM)----
fig4_metadata <- metadata %>% 
  filter(group == "C" & exp_num %in% c("M7","M9")| #Only use C mice from these experiments. Allocated groups to figures based on paper outline.
         group == "CWM" & exp_num %in% c("M6","M7", "M9")|
         group == "RM" & exp_num %in% c("M7","M9")|
         group == "FRM" & exp_num %in% c("M9")) %>% 
  mutate(group=factor(group, levels=c("C", "CWM", "RM", "FRM")))  # Make sure group is treated as a factor


# of mice represented in the figure
fig4_mice <- length(unique(fig4_metadata$m_id_unique)) 
# 48 mice total for figure 1

#Note baseline weight for each group of mice (based on the earliest timepoint recorded for each experiment)----
baseline <- fig4_metadata %>% #Baseline weight was taken at day -5 for groups C, WM, and WMC
  filter(group == "C" & day == -2| #12 mice in C group
         group == "CWM" & day == -15| #6 mice in CWM group with baseline at day -15
         group == "CWM" & day == -2 & exp_num %in% c("M7", "M9")| #12 mice in CWM group with baseline at day -2
         group == "RM" & day == -2| #12 mice in RM group
         group == "FRM" & day == -2) %>% #6 mice in FRM group
  mutate(baseline_weight = weight) %>% #This column represents the initial weight that was recorded for each mouse
  select(m_id_unique, baseline_weight) #Will use m_id_unique to join baseline_weights to fig4_metadata

#Make a new column that represents weight_change from baseline_weight----
fig4_weightdata <- inner_join(fig4_metadata, baseline, by = "m_id_unique") %>% #Join baseline weight to fig4_metadata
  group_by(m_id_unique, day) %>% #Group by each unique mouse and experiment day
  mutate(weight_change = weight-baseline_weight) %>% #Make a new column that represents the change in weight from baseline (all weights recorded in grams)
  ungroup() %>% 
  filter(!is.na(weight)) #drop rows with NA values for fig4_weightdata. 744 samples including NAs, 744 samples after excluding NAs

#C. diff CFU plot----
#Narrow fig4_metadata to just timepoints relevant to C. difficile CFU tracking (Anything on or after day 0)
fig4_cfudata <- fig4_metadata %>% 
  filter(day > -1)
fig4_cfu_na <- sum(is.na(fig4_cfudata$avg_cfu)) #14 samples with NA values. Represent times when we either did not collect stool samples or weren't able to get a stool sample from a particular mouse
#Drop rows with NA values for fig4_cfu:
fig4_cfudata <- fig4_cfudata %>% 
  filter(!is.na(avg_cfu))

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
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                       limits = c(-1, 31)) +
    geom_hline(yintercept = 100, linetype=2) + #Line that represents our limit of detection when quantifying C. difficile CFU by plating
    geom_text(x = 11, y = 104, color = "black", label = "LOD") + #Label for line that represents our limit of detection when quantifying C. difficile CFU by plating
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9))+ #Scientific notation labels for y-axis
    theme_classic()
}
fig4_cfu <- plot_cfu(fig4_cfudata)
save_plot(filename = "results/figures/fig4_cfu.png", fig4_cfu, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)


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
    scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20, 25, 30),
                       limits = c(-16, 31)) +
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
    scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20, 25, 30),
                       limits = c(-16, 31)) +
    theme_classic()
}

#Function to plot weight that only plots the mean of each group from D0 through D10 and no points for individual mice. Argument = dataframe you want to plot.
plot_weight_10d <- function(df){
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
    scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8, 10),
                       limits = c(-3, 11)) +
    theme_classic()
}

fig4_weight <- plot_weight(fig4_weightdata) 
save_plot(filename = "results/figures/fig4_weight.png", fig4_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
fig4v2_weight <- plot_weight_simple(fig4_weightdata)
save_plot(filename = "results/figures/fig4v2_weight.png", fig4v2_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
fig4_weight_10d <- plot_weight_10d(fig4_weightdata)
save_plot(filename = "results/figures/fig4_weight_10d.png", fig4_weight_10d, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
kruskal_wallis_cfu <- fig4_cfudata %>% 
  filter(day %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15)) %>%  #only test days that we have CFU data for #Only have cfu for CWM group on D20,25,30, exclude those days
  group_by(day) %>% 
  do(tidy(kruskal.test(avg_cfu~factor(group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) 
#Timepoints where C. diff CFU is significantly different across the sources of mice
sig_C.diff_CFU_timepoints <- kruskal_wallis_cfu %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day)
