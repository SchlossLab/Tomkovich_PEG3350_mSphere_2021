source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Define color scheme for this figure----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMR")
color_labels <- c( "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 + 10-day recovery")

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
sig_cfu_days <- cfu_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  pull(day)

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
sig_weight_days <- weight_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  pull(day)

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

#Dataframe of cfu data for just the initial 10 days of the experiment
cfudata_10dsubset <- cfudata %>%
  filter(day < 12) #only include data through day 10

#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- cfu_kruskal_wallis_adjust %>%
  filter(day < 12) %>% #Only include results through day 10 for this plot
  filter(p.value.adj <= 0.05) %>%
  pull(day)
y_position <- max(cfudata$avg_cfu)
label <- cfu_kruskal_wallis_adjust %>%
  filter(day < 12) %>% #Only include results through day 10 for this plot
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>%
  pull(p.signif)

#Plot cfu for just the inital 10days
cfu_10d <- plot_cfu_data(cfudata_10dsubset) +
      scale_x_continuous(breaks = c(0:10),
                         limits = c(-1, 11),
                         minor_breaks = c(-.5:10.5))
save_plot(filename = "results/figures/5_days_PEG_cfu_10d.png", cfu_10d, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Statistical annotation labels based on adjusted kruskal-wallis p-values for all timepoints:
x_annotation <- sig_cfu_days
y_position <- max(cfudata$avg_cfu)
label <- cfu_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>%
  pull(p.signif)

#Plot of cfu data for all days of the experiment
cfu <- plot_cfu_data(cfudata) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),
                     limits = c(-1, 31),
                     minor_breaks = c(-.5:10.5, 11.5, 12.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5)) #only show grey lines separating days on days with statistically sig points
save_plot(filename = "results/figures/5_days_PEG_cfu.png", cfu, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Weight change plot----

#Dataframe of weight data for days -15 through 10 of the experiment:
weight_subset <- weightdata %>%
  filter(day < 12)

#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- weight_kruskal_wallis_adjust %>%
  filter(day < 12) %>% #Only include results through day 10 for this plot
  filter(p.value.adj <= 0.05) %>%
  pull(day)
x_annotation == sig_weight_days #All the days where weight change significantly varied across groups of mice occured within the first 10 days post-infection
y_position <- max(weightdata$weight_change)
label <- weight_kruskal_wallis_adjust %>%
  filter(day < 12) %>% #Only include results through day 10 for this plot
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>%
  pull(p.signif)

#Plot of weight data for days -15 through 10 of the experiment:
weight_subset_plot <- plot_weight(weight_subset) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10),
                     limits = c(-16, 11),
                     minor_breaks = c(-15.5:10.5)) 
save_plot(filename = "results/figures/5_days_PEG_weight_subset.png", weight_subset_plot, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)


#Plot of weight data for all days of the experiment:
#Note don't need to redo statistical annotations since there were no significant differences past 1 day post-infection
weight_plot <- plot_weight(weightdata) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15, 20, 25, 30),
                     limits = c(-16, 31),
                     minor_breaks = c(-15.5:10.5, 11.5, 12.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5)) #only show grey lines around days on days with points)

#Plots with just the median lines for each group
v2_weight_subset <- plot_weight_medians(weight_subset) +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10),
                     limits = c(-16, 11),
                     minor_breaks = c(-15.5:10.5)) #only show grey lines separating days on days with statistically sig points)
save_plot(filename = "results/figures/5_days_PEGv2_weight_subset.png", v2_weight_subset, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

