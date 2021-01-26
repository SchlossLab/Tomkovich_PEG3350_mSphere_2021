source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Define color scheme for this figure----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "CWM", "FRM", "RM")
color_labels <- c( "Clind.", "Clind. + 1-day PEG 3350", "Clind. + 3-day recovery + 1-day PEG 3350 + FMT", "Clind. + 3-day recovery + 1-day PEG 3350")

#Subset metadata to relevant groups and experiments (C, CWM, RM, FRM)----
metadata <- post_cdi_PEG_subset(metadata) %>% 
  filter(!sample_type %in% c("cecum", "distal_colon", "proximal_colon")) %>% #Get rid of rows corresponding to tissue samples in the metadata as these will create duplicate values for mice at timepoints where tissues were also collected
  mutate(day = as.integer(day))  #Day variable transformed to integer 

# of mice represented in the figure
mice <- length(unique(metadata$unique_mouse_id))
# 48 mice total for 5_days_PEG figure

sum <- metadata %>%
  group_by(group) %>%
  count(day)

#C. difficile CFU dataframe----
#Narrow metadata to just timepoints relevant to C. difficile CFU tracking (Anything on or after day 0)
cfudata <- metadata %>%
  filter(day > -1)
cfu_na <- sum(is.na(cfudata$avg_cfu)) #14 samples with NA values. Represent times when we either did not collect stool samples or weren't able to get a stool sample from a particular mouse
#Drop rows with NA values for cfu:
cfudata <- cfudata %>%
  filter(!is.na(avg_cfu))

#Weight change dataframe----
#Note baseline weight for each group of mice (based on the earliest timepoint recorded for each experiment)----
baseline <- metadata %>% #Baseline weight was taken at day -5 for groups C, WM, and WMC
  filter(group == "C" & day == -2| #12 mice in C group
         group == "CWM" & day == -15| #6 mice in CWM group with baseline at day -15
         group == "CWM" & day == -2 & exp_num %in% c("M7", "M9")| #12 mice in CWM group with baseline at day -2
         group == "RM" & day == -2| #12 mice in RM group
         group == "FRM" & day == -2) %>% #6 mice in FRM group
  mutate(baseline_weight = weight) %>% #This column represents the initial weight that was recorded for each mouse
  select(unique_mouse_id, baseline_weight) #Will use unique_mouse_id to join baseline_weights to metadata

#Make a new column that represents weight_change from baseline_weight----
weightdata <- inner_join(metadata, baseline, by = "unique_mouse_id") %>% #Join baseline weight to metadata
  group_by(unique_mouse_id, day) %>% #Group by each unique mouse and experiment day
  mutate(weight_change = weight-baseline_weight) %>% #Make a new column that represents the change in weight from baseline (all weights recorded in grams)
  ungroup() %>%
  filter(!is.na(weight)) #drop rows with NA values for weightdata. 744 samples including NAs, 744 samples after excluding NAs

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Statiscal analysis of C. difficile CFU data----
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
cfu_kruskal_wallis <- cfudata %>%
  filter(day %in% c(0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 15)) %>%  #only test days that we have CFU data for at least 3 groups
  select(day, group, avg_cfu) %>%
  group_by(day) %>%
  nest() %>%
  mutate(model=map(data, ~kruskal.test(x=.x$avg_cfu, g=as.factor(.x$group)) %>% tidy())) %>%
  mutate(median = map(data, get_cfu_median)) %>%
  unnest(c(model, median)) %>%
  ungroup()
#Adjust p-values for testing multiple days and write results to table:
cfu_kruskal_wallis_adjust <- cfu_kruskal_wallis %>%
  select(day, statistic, p.value, parameter, method, C, CWM, RM, FRM) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv("data/process/post_CDI_PEG_cfu_stats_all_days.tsv")

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
  write_tsv("data/process/post_CDI_PEG_cfu_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
cfu_plot_format_stats <- cfu_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method,-C, -CWM, -RM, -FRM) %>%
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>%
  bind_rows()

#Statistical analysis of mouse weight change data----
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
weight_kruskal_wallis <- weightdata %>%
  filter(day %in% c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15)) %>%  #only test days that we have weight data for at least 3 groups
  select(day, group, weight_change) %>%
  group_by(day) %>%
  nest() %>%
  mutate(model=map(data, ~kruskal.test(x=.x$weight_change, g=as.factor(.x$group)) %>% tidy())) %>%
  mutate(median = map(data, get_weight_median)) %>%
  unnest(c(model, median)) %>%
  ungroup()
#Adjust p-values for testing multiple days and write results to table:
weight_kruskal_wallis_adjust <- weight_kruskal_wallis %>%
  select(day, statistic, p.value, parameter, method, C, CWM, RM, FRM) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj) %>%
  write_tsv("data/process/post_CDI_PEG_weight_stats_all_days.tsv")

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
  write_tsv("data/process/post_CDI_PEG_weight_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
weight_plot_format_stats <- weight_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method,-C, -CWM, -RM, -FRM) %>%
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>%
  bind_rows()

#Plots of CFU and weight data----

#C. diff CFU plot
#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- cfu_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  pull(day)
y_position <- max(cfudata$avg_cfu) + 500000000
label <- kw_label(cfu_kruskal_wallis_adjust)

cfu <- plot_cfu_data(cfudata) +
  scale_x_continuous(breaks = c(0:10, 15, 20, 25, 30),
                     limits = c(-1, 31),
                     minor_breaks = c(-.5:10.5, 14.5, 15.5, 19.5, 20.5, 24.5, 25.5, 29.5, 30.5))+
  theme(legend.position = "none")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = Inf, fill = "#88419d", alpha = .15)+ #shade to indicate PEG treatment in Clind + 1-day PEG group
  annotate("rect", xmin = 3, xmax = 4, ymin = 0, ymax = Inf, fill = "#225ea8", alpha = .15)+ #shade to inidcate PEG treatment in Clind + 3-day recovery + 1-day PEG + FMT/PBS
  geom_vline(aes(xintercept = 3, alpha = .4), color = "#f768a1")+ #Lines to indicate the groups with FMTs and PBS both received PEG during this period
  geom_vline(aes(xintercept = 4, alpha = .4), color = "#f768a1") +
  geom_rect(aes(xmin = 15.75, xmax = 31, ymin = 104, max = 2200), fill = "white")+ #blank background for legend
  geom_rect(aes(xmin = 16, xmax = 17.5, ymin = 900, max = 1800), fill = "#88419d", alpha = .2)+ #box for Clind + 1-day PEG group legend
  geom_rect(aes(xmin = 16, xmax = 17.5, ymin = 200, max = 400), fill = "#225ea8", color = "#f768a1", alpha = .2)+ #box for Clind + 3-day recovery + 1-day PEG + FMT/PBS group legend
  geom_text(aes(x= 24, y = 1300, label = "PEG Treatment for Clind. + 1-day PEG 3350 Group"), size = 2.25)+
  geom_text(aes(x= 25, y = 300, label = " PEG Treatment for Clind. + 3-day recovery + 1-day PEG 3350 (+ FMT) Groups"), size = 2.25)
save_plot(filename = "results/figures/post_CDI_PEG_cfu.png", cfu, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Plot of just a subset of data (through day 15)
cfu_subset <- plot_cfu_data(cfudata %>% filter(day < 16)) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 15),
                     limits = c(-1, 16),
                     minor_breaks = c(-.5:10.5, 14.5, 15.5)) #only show grey lines indicating day for stat sig days
save_plot(filename = "results/figures/post_CDI_PEG_cfu_subset.png", cfu_subset, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Weight change plot----
x_annotation <- weight_kruskal_wallis_adjust %>%
  filter(p.value.adj <= 0.05) %>%
  pull(day)
y_position <- max(weightdata$weight_change)
label <- kw_label(weight_kruskal_wallis_adjust)


weight <- plot_weight(weightdata %>% filter(day %in% c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15))) + #Narrow weight data to just timepoints where we have data for at least 3 groups
  scale_x_continuous(breaks = c(-2, -4, -2, 0, 2, 4, 6, 8, 10, 15),
                     limits = c(-3, 16),
                     minor_breaks = c(-2.5:10.5, 14.5, 15.5)) 
save_plot(filename = "results/figures/post_CDI_PEG_weight.png", weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#Show just median lines for each group
y_position <- 2 #Change for the plot showing just the median lines
v2_weight <- plot_weight_medians(weightdata %>% filter(day %in% c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15))) + #Narrow weight data to just timepoints where we have data for at least 3 groups
  scale_x_continuous(breaks = c(-2, -4, -2, 0, 2, 4, 6, 8, 10, 15),
                     limits = c(-3, 16),
                     minor_breaks = c(.5:8.5)) #only show grey lines indicating day for stat sig days
save_plot(filename = "results/figures/post_CDI_PEGv2_weight.png", v2_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

weight_10d <- plot_weight_medians(weightdata %>% filter(day %in% c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))) + #Narrow weight data to just timepoints where we have data for at least 3 groups
  scale_x_continuous(breaks = c(-2, -4, -2, 0, 2, 4, 6, 8, 10),
                     limits = c(-3, 11),
                     minor_breaks = c(.5:8.5)) #only show grey lines indicating day for stat sig days
save_plot(filename = "results/figures/post_CDI_PEG_weight_10d.png", weight_10d, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
