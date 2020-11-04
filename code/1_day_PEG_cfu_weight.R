source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Define color scheme for this figure----
color_scheme <- c("#225ea8", "#238b45", "#88419d") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c( "1RM1", "C", "M1")
color_labels <- c( "1-day PEG 3350 + 1-day recovery", "Clind.", "1-day PEG 3350")

#Subset metadata to relevant groups and experiments (C, 1RM1, M1) for this subset of mice----
metadata <- one_day_PEG_subset(metadata)


# of mice represented in the figure
mice <- length(unique(metadata$unique_mouse_id))
# 18 mice total for figure 3

#C. difficile CFU dataframe----
#Narrow metadata to just timepoints relevant to C. difficile CFU tracking (Anything on or after day 0)
cfudata <- metadata %>%
  filter(day > -1)
cfu_na <- sum(is.na(cfudata$avg_cfu)) #53 samples with NA values. 5 samples, where we weren't able to get a stool sample. Rest of NAs are from timepoints after D15 which is the day we stopped tracking C. diff CFU
#Drop rows with NA values for cfu:
cfudata <- cfudata %>%
  filter(!is.na(avg_cfu)) #181 samples total

#Weight change dataframe----
#Note baseline weight for each group of mice (based on the earliest timepoint recorded for each experiment)----
baseline <- metadata %>% #Baseline weight was taken at day -5 for groups C, WM, and WMC
  filter(group == "C" & day == -15| #6 mice in C group
           group == "1RM1" & day == -2| #6 mice in 1RM1 group
           group == "M1" & day == -11) %>%  #6 mice in M1 group
  mutate(baseline_weight = weight) %>% #This column represents the initial weight that was recorded for each mouse
  select(unique_mouse_id, baseline_weight) #Will use unique_mouse_id to join baseline_weights to metadata

#Make a new column that represents weight_change from baseline_weight----
weightdata <- inner_join(metadata, baseline, by = "unique_mouse_id") %>% #Join baseline weight to metadata
  group_by(unique_mouse_id, day) %>% #Group by each unique mouse and experiment day
  mutate(weight_change = weight-baseline_weight) %>% #Make a new column that represents the change in weight from baseline (all weights recorded in grams)
  ungroup() %>%
  filter(!is.na(weight)) #drop rows with NA values for weightdata. 378 samples including NAs, 306 samples after excluding NAs

#Statistical analysis of C. difficile CFU data----
set.seed(19760620) #Same seed used for mothur analysis
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
cfu_kruskal_wallis <- cfudata %>%
  filter(day %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8)) %>%  #only test days that we have cfu data for at least 3 groups
  select(day, group, avg_cfu) %>%
  group_by(day) %>%
  nest() %>%
  mutate(model=map(data, ~kruskal.test(x=.x$avg_cfu, g=as.factor(.x$group)) %>% tidy())) %>%
  mutate(median = map(data, get_cfu_median)) %>%
  unnest(c(model, median)) %>%
  ungroup()

#Adjust p-values for testing multiple days and write results to table:
cfu_kruskal_wallis_adjust <- cfu_kruskal_wallis %>%
  select(day, statistic, p.value, parameter, method, `1RM1`, C, M1) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv("data/process/1_day_PEG_cfu_stats_all_days.tsv")

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
  write_tsv("data/process/1_day_PEG_cfu_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
cfu_plot_format_stats <- cfu_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method, -C, -M1, -`1RM1`) %>%
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>%
  bind_rows()

#Statistical analysis of mouse weight change data----
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
weight_kruskal_wallis <- weightdata %>%
  filter(day %in% c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) %>%  #only test days that we have weight data for at least 3 groups
  select(day, group, weight_change) %>%
  group_by(day) %>%
  nest() %>%
  mutate(model=map(data, ~kruskal.test(x=.x$weight_change, g=as.factor(.x$group)) %>% tidy())) %>%
  mutate(median = map(data, get_weight_median)) %>%
  unnest(c(model, median)) %>%
  ungroup()

#Adjust p-values for testing multiple days and write results to table:
weight_kruskal_wallis_adjust <- weight_kruskal_wallis %>%
  select(day, statistic, p.value, parameter, method, -C, -M1, -`1RM1`) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv("data/process/1_day_PEG_weight_stats_all_days.tsv")

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
  write_tsv("data/process/1_day_PEG_weight_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
weight_plot_format_stats <- weight_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method,-C, -M1, -`1RM1`) %>%
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>%
  bind_rows()

#Plots of CFU and weight data----

#Plot of CFU data----
#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- cfu_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  pull(day)
y_position <- max(cfudata$avg_cfu) + 100000000
label <- cfu_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>%
  pull(p.signif)

cfu <- plot_cfu_data(cfudata) +
  scale_x_continuous(breaks = c(0:10),
                     limits = c(-1, 11),
                     minor_breaks = c(-.5:10.5))

save_plot(filename = "results/figures/1_day_PEG_cfu.png", cfu, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Plot of weight change data----
#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- weight_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  pull(day)
y_position <- max(weightdata$weight_change)
label <- weight_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>%
  pull(p.signif)

#Narrow data to just the timepoints tested in statistical analysis, when we have weight data for all 3 groups
weightdata_subset <- weightdata %>%
  filter(day %in% c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

weight <- plot_weight(weightdata_subset) +
  scale_x_continuous(breaks = c(-2:10),
                     limits = c(-2.5, 10.5),
                     minor_breaks = c(-2.5:10.5))
save_plot(filename = "results/figures/1_day_PEG_weight.png", weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

y_position <- 2 #Change for the plot showing just the median lines
v2_weight <- plot_weight_medians(weightdata_subset)+
  scale_x_continuous(breaks = c(-2:10),
                     limits = c(-3, 11),
                     minor_breaks = c(-2.5:10.5))
save_plot(filename = "results/figures/1_day_PEGv2_weight.png", v2_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

