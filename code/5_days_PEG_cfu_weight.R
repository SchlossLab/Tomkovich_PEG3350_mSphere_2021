source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Define color scheme for this figure----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMR")
color_labels <- c( "Clind.", "5-day PEG", "5-day PEG + Clind.", "5-day PEG + 10-day recovery")

#Subset metadata to relevant groups and experiments (WM, WMC, WMR, C)----
metadata <- five_day_PEG_subset(metadata) %>% 
  filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon"))  #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected

# of mice represented in the figure
mice <- length(unique(metadata$unique_mouse_id))
# 62 mice total for 5_days_PEG figure

#C. difficile CFU dataframe----
#Narrow metadata to just timepoints relevant to C. difficile CFU tracking (Anything on or after day 0)
cfudata <- metadata %>%
  filter(day > -1)
cfu_na <- sum(is.na(cfudata$avg_cfu)) #182 samples with NA values. Represent times when we either did not collect stool samples, weren't able to get a stool sample from a particular mouse, weren't able to plate the sample we did collect immediately after due to chamber issues or time constraints, or the mouse died early
#Drop rows with NA values for cfu:
cfudata <- cfudata %>%
  filter(!is.na(avg_cfu))

#Weight change dataframe----
#Note baseline weight for each group of mice (based on the earliest timepoint recorded for each experiment)----
baseline <- metadata %>% #Baseline weight was taken at day -5 for groups C, WM, and WMC
  filter(group == "C" & day == -5| #20 mice in C group
         group == "WM" & day == -5| #21 mice in WM group
         group == "WMC" & day == -5| #9 mice in WMC group
         group == "WMR" & day == -15) %>% #12 mice in WMR group, baseline weight was taken at day -15
  mutate(baseline_weight = weight) %>% #This column represents the initial weight that was recorded for each mouse
  select(unique_mouse_id, baseline_weight) #Will use unique_mouse_id to join baseline_weights to metadata

#Make a new column that represents weight_change from baseline_weight
weightdata <- inner_join(metadata, baseline, by = "unique_mouse_id") %>% #Join baseline weight to metadata
  group_by(unique_mouse_id, day) %>% #Group by each unique mouse and experiment day
  mutate(weight_change = weight-baseline_weight) %>% #Make a new column that represents the change in weight from baseline (all weights recorded in grams)
  ungroup() %>%
  filter(!is.na(weight)) #drop rows with NA values for weightdata. 1040 samples including NAs, 870 samples after excluding NAs

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis
#Shapiro-Wilk test to see if cfu and weight change data is normally distributed:
#Note: p-value > 0.05 means the data is normally distributed
shapiro.test(cfudata$avg_cfu) #p-value < 2.2e-16
shapiro.test(weightdata$weight_change) #p-value = 1.485e-09
#Since p-value < 0.05 for both variables, we will use non-parametric tests

#Statiscal analysis of C. difficile CFU data----
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
cfu_kruskal_wallis <- cfudata %>%
  filter(day %in% c(0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 15, 20, 25, 30)) %>%  #only test days that we have CFU data for #Only have cfu for WMR group on D7, exclude that day
  select(day, group, avg_cfu) %>%
  group_by(day) %>%
  nest() %>%
  mutate(model=map(data, ~kruskal.test(x=.x$avg_cfu, g=as.factor(.x$group)) %>% tidy())) %>%
  mutate(median = map(data, get_cfu_median)) %>%
  unnest(c(model, median)) %>%
  ungroup()
#Adjust p-values for testing multiple days and write results to table:
cfu_kruskal_wallis_adjust <- cfu_kruskal_wallis %>%
  select(day, statistic, p.value, parameter, method, C, WM, WMC, WMR) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv("data/process/5_days_PEG_cfu_stats_all_days.tsv")

#Timepoints where C. difficile CFU is significantly different across the groups of mice after BH adjustment of p-values:
sig_cfu_days <- pull_sig_days(cfu_kruskal_wallis_adjust)

#Perform pairwise Wilcoxan rank sum tests for days that were significant by Kruskal-Wallis test
cfu_stats_pairwise <- cfu_kruskal_wallis %>%
  filter(day %in% sig_cfu_days) %>% #only perform pairwise tests for days that were significant
  group_by(day) %>%
  mutate(model=map(data, ~pairwise.wilcox.test(x=.x$avg_cfu, g=as.factor(.x$group), p.adjust.method="BH") %>%
                     tidy() %>%
                     mutate(compare=paste(group1, group2, sep="-")) %>%
                     select(-group1, -group2) %>%
                     pivot_wider(names_from=compare, values_from=p.value)
  )
  ) %>%
  unnest(model) %>%
  select(-data, -parameter, -statistic) %>%
  write_tsv("data/process/5_days_PEG_cfu_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
cfu_plot_format_stats <- cfu_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method,-C, -WM, -WMC, -WMR) %>%
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>%
  bind_rows()

#Statistical analysis of mouse weight change data----
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
weight_kruskal_wallis <- weightdata %>%
  filter(day %in% c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30)) %>%  #only test days that we have weight data for at least 3 groups
  select(day, group, weight_change) %>%
  group_by(day) %>%
  nest() %>%
  mutate(model=map(data, ~kruskal.test(x=.x$weight_change, g=as.factor(.x$group)) %>% tidy())) %>%
  mutate(median = map(data, get_weight_median)) %>%
  unnest(c(model, median)) %>%
  ungroup()
#Adjust p-values for testing multiple days and write results to table:
weight_kruskal_wallis_adjust <- weight_kruskal_wallis %>%
  select(day, statistic, p.value, parameter, method, C, WM, WMC, WMR) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv("data/process/5_days_PEG_weight_stats_all_days.tsv")

#Timepoints where C. difficile CFU is significantly different across the groups of mice after BH adjustment of p-values:
sig_weight_days <- pull_sig_days(weight_kruskal_wallis_adjust)

#Perform pairwise Wilcoxan rank sum tests for days that were significant by Kruskal-Wallis test
weight_stats_pairwise <- weight_kruskal_wallis %>%
  filter(day %in% sig_weight_days) %>% #only perform pairwise tests for days that were significant
  group_by(day) %>%
  mutate(model=map(data, ~pairwise.wilcox.test(x=.x$weight_change, g=as.factor(.x$group), p.adjust.method="BH") %>%
                     tidy() %>%
                     mutate(compare=paste(group1, group2, sep="-")) %>%
                     select(-group1, -group2) %>%
                     pivot_wider(names_from=compare, values_from=p.value)
  )
  ) %>%
  unnest(model) %>%
  select(-data, -parameter, -statistic) %>%
  write_tsv("data/process/5_days_PEG_weight_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
weight_plot_format_stats <- weight_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method,-C, -WM, -WMC, -WMR) %>%
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>%
  bind_rows()

#Plots of CFU and weight data----
#Transform day column variable from character to integer variable
cfudata <- cfudata %>% 
  mutate(day = as.integer(day))

#Dataframe of cfu data for just the initial 10 days of the experiment
cfudata_10dsubset <- cfudata %>%
  mutate(day = as.integer(day)) %>% #transform day to integer variable
  filter(day < 12) #only include data through day 10

cfu_kruskal_wallis_adjust <- cfu_kruskal_wallis_adjust %>% 
  mutate(day = as.integer(day)) #transform day to integer variable
  
#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- cfu_kruskal_wallis_adjust %>%
  filter(day < 12) %>% #Only include results through day 10 for this plot
  filter(p.value.adj <= 0.05) %>%
  pull(day)
y_position <- max(cfudata$avg_cfu)
label <- kw_label(cfu_kruskal_wallis_adjust %>% filter(day < 12)) 
 #Only include results through day 10 for this plot

#Plot cfu for just the inital 10days
cfu_10d <- plot_cfu_data(cfudata_10dsubset %>% 
                           filter(day %in% c("-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) + #Filter to just include timepoints that will be plotted
      scale_x_continuous(breaks = c(0:10),
                         limits = c(-1, 11),
                         minor_breaks = c(-.5:10.5))
save_plot(filename = "results/figures/5_days_PEG_cfu_10d.png", cfu_10d, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Statistical annotation labels based on adjusted kruskal-wallis p-values for all timepoints:
x_annotation <- cfu_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  pull(day)
y_position <- max(cfudata$avg_cfu)
label <- kw_label(cfu_kruskal_wallis_adjust)

#Plot of cfu data for all days of the experiment
cfu <- plot_cfu_data(cfudata %>% 
                       filter(day %in% c("-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15", "20", "25", "30"))) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                     limits = c(-1, 31),
                     minor_breaks = c(-.5:10.5, 11.5, 12.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5))+ #only show grey lines separating days on days with statistically sig points
  theme(legend.position = "none")
save_plot(filename = "results/figures/5_days_PEG_cfu.png", cfu, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Plot of just the cfu data for the WMR group (5-day PEG + 10-day recovery)
wmr_cfu <- cfudata %>% 
  filter(group == "WMR") %>% 
  filter(day %in% c("-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15", "20", "25", "30"))
wmr_median_summary <- wmr_cfu %>%
  group_by(group, day) %>%
  summarize(median_avg_cfu = median(avg_cfu, na.rm = TRUE))
wmr_cfu_plot <-  ggplot(NULL) +
    geom_point(wmr_cfu, mapping = aes(x = day, y = avg_cfu, color= group, fill = group), alpha = 0.7, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(wmr_median_summary, mapping = aes(x = day, y = median_avg_cfu, group = group, color = group), alpha = 0.6, size = 1.5) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    labs(x = "Days post-challenge", y = "CFU/g feces") +
    scale_y_log10(breaks = c(100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10), 
                  labels = c('10^2', '10^3', '10^4', '10^5', '10^6', '10^7', '10^8', '10^9', '10^10')) + # scale y axis log10 and label 10^x
    geom_hline(yintercept = 100, linetype=2) + #Line that represents our limit of detection when quantifying C. difficile CFU by plating
    geom_text(x = 11, y = 104, color = "black", label = "LOD") + #Label for line that represents our limit of detection when quantifying C. difficile CFU by plating
    theme(text = element_text(size = 16))+  # Change font size for entire plot
    theme_classic()+
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                     limits = c(-1, 31),
                     minor_breaks = c(-.5:10.5, 11.5, 12.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5))+ #only show grey lines separating days on days with statistically sig points
    theme(legend.position = "none",
          axis.text.y = element_markdown(size = 12),
          legend.key= element_rect(colour = "transparent", fill = "transparent"),
          text = element_text(size = 16), # Change font size for entire plot
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days)
save_plot(filename = "results/figures/5_days_PEG_cfu_WMR.png", wmr_cfu_plot, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Make a second version where each unique mouse = a different colored line
#Plot of just the cfu data for the WMR group (5-day PEG + 10-day recovery)
wmr_cfu_mice <- cfudata %>% 
  filter(group == "WMR") %>% 
  filter(day %in% c("-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15", "20", "25", "30")) %>% 
  filter(!duplicated(unique_mouse_id)) %>% #Remove duplicate mouse ids
  mutate(unique_mouse_id = factor(unique_mouse_id, levels = unique(as.factor(unique_mouse_id)))) %>% 
  pull(unique_mouse_id)
color_mice <- c("5_M5",  "6_M5",  "7_M5",  "8_M5",  "9_M5",  "10_M5", "7_M6",  "9_M6",  "10_M6", "11_M6", "12_M6")
color_mice_values <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
                       "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#9B870C")
color_mice_labels <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10 , 11)

wmr_cfu_plot_indiv <-  wmr_cfu %>% 
  ggplot()+
  geom_line(mapping = aes(x = day, y = avg_cfu, group = unique_mouse_id, color = unique_mouse_id), alpha = 0.8, size = 1.5) +
  scale_colour_manual(name="Mouse",
                      values=color_mice_values,
                      breaks=color_mice,
                      labels=color_mice_labels)+
  labs(x = "Days post-challenge", y = "CFU/g feces") +
  scale_y_log10(breaks = c(100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10), 
                labels = c('10^2', '10^3', '10^4', '10^5', '10^6', '10^7', '10^8', '10^9', '10^10')) + # scale y axis log10 and label 10^x
  geom_hline(yintercept = 100, linetype=2) + #Line that represents our limit of detection when quantifying C. difficile CFU by plating
  geom_text(x = 11, y = 104, color = "black", label = "LOD") + #Label for line that represents our limit of detection when quantifying C. difficile CFU by plating
  theme(text = element_text(size = 16))+  # Change font size for entire plot
  theme_classic()+
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                     limits = c(-1, 31),
                     minor_breaks = c(-.5:10.5, 11.5, 12.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5))+ #only show grey lines separating days on days with statistically sig points
  guides(colour = guide_legend(nrow = 1))+#Limit number of rows in the legend
  theme(legend.position = "bottom",
        axis.text.y = element_markdown(size = 12),
        text = element_text(size = 16), # Change font size for entire plot
        axis.ticks.x = element_blank(),
        legend.key= element_rect(colour = "transparent", fill = "transparent"),
        panel.grid.minor.x = element_line(size = 0.4, color = "grey"))#Add gray lines to clearly separate symbols by days)
save_plot(filename = "results/figures/5_days_PEG_cfu_WMR_indiv.png", wmr_cfu_plot_indiv, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Weight change plot----

#Dataframe of weight data for days -15 through 10 of the experiment:
weight_subset <- weightdata %>%
  filter(day < 12)

weight_kruskal_wallis_adjust <- weight_kruskal_wallis_adjust %>% 
  mutate(day = as.integer(day))

#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- weight_kruskal_wallis_adjust %>%
  filter(day < 12) %>% #Only include results through day 10 for this plot
  filter(p.value.adj <= 0.05) %>%
  pull(day)
x_annotation == sig_weight_days #All the days where weight change significantly varied across groups of mice occured within the first 10 days post-infection
y_position <- max(weightdata$weight_change)
label <- kw_label(weight_kruskal_wallis_adjust %>%
                    filter(day < 12)) #Only include results through day 10 for this plot

weight_subset <- weight_subset %>% 
  mutate(day = as.integer(day))
weightdata <- weightdata %>% 
  mutate(day = as.integer(day))

#Plot of weight data for days -15 through 10 of the experiment:
weight_subset_plot <- plot_weight(weight_subset) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10),
                     limits = c(-16, 11),
                     minor_breaks = c(-15.5:10.5)) 
save_plot(filename = "results/figures/5_days_PEG_weight_subset.png", weight_subset_plot, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Plots with just the median lines for each group
v2_weight_subset <- plot_weight_medians(weight_subset) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10),
                     limits = c(-16, 11),
                     minor_breaks = c(-15.5:10.5)) #only show grey lines separating days on days with statistically sig points)
save_plot(filename = "results/figures/5_days_PEGv2_weight_subset.png", v2_weight_subset, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Plot of weight data for all days of the experiment:
#Note don't need to redo statistical annotations since there were no significant differences past 10 day post-infection
#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- weight_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  pull(day)
x_annotation == sig_weight_days #All the days where weight change significantly varied across groups of mice occured within the first 10 days post-infection
y_position <- max(weightdata$weight_change)
label <- kw_label(weight_kruskal_wallis_adjust) 

weight_plot <- plot_weight(weightdata) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20, 25, 30),
                     limits = c(-16, 31),
                     minor_breaks = c(-15.5:10.5, 11.5, 12.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5)) #only show grey lines around days on days with points)
save_plot(filename = "results/figures/5_days_PEG_weight.png", weight_plot , base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Plots with just the median lines for each group
v2_weight_plot <- plot_weight_medians(weightdata) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20, 25, 30),
                     limits = c(-16, 31),
                     minor_breaks = c(-15.5:10.5, 11.5, 12.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5)) #only show grey lines around days on days with points)
save_plot(filename = "results/figures/5_days_PEG_weight_median.png", v2_weight_plot, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
