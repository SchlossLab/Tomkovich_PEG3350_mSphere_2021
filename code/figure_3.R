source("code/functions.R") #Reads in metadata

#Define color scheme for this figure----
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

#C. difficile CFU dataframe----
#Narrow fig3_metadata to just timepoints relevant to C. difficile CFU tracking (Anything on or after day 0)
fig3_cfudata <- fig3_metadata %>% 
  filter(day > -1)
fig3_cfu_na <- sum(is.na(fig3_cfudata$avg_cfu)) #53 samples with NA values. 5 samples, where we weren't able to get a stool sample. Rest of NAs are from timepoints after D15 which is the day we stopped tracking C. diff CFU
#Drop rows with NA values for fig3_cfu:
fig3_cfudata <- fig3_cfudata %>% 
  filter(!is.na(avg_cfu)) #181 samples total

#Weight change dataframe----
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

#Statiscal analysis of C. difficile CFU data----
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
cfu_kruskal_wallis <- fig3_cfudata %>% 
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
  write_tsv("data/process/fig3_cfu_stats_all_days.tsv")

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
  write_tsv("data/process/fig3_cfu_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
cfu_plot_format_stats <- cfu_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method, -C, -M1, -`1RM1`) %>% 
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>% 
  bind_rows()

#Statistical analysis of mouse weight change data----
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
weight_kruskal_wallis <- fig3_weightdata %>% 
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
  write_tsv("data/process/fig3_weight_stats_all_days.tsv")

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
  write_tsv("data/process/fig3_weight_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
weight_plot_format_stats <- weight_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method,-C, -M1, -`1RM1`) %>% 
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>% 
  bind_rows()

#Plots of CFU and weight data----
#Function to plot cfu data
#Arguments: df = dataframe to plot
#When using the function can add ggplot lines to specify x axis scale
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

#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- cfu_kruskal_wallis_adjust %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day)
y_position <- max(fig3_cfudata$avg_cfu) + 100000000
label <- cfu_kruskal_wallis_adjust %>% 
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>% 
  pull(p.signif)

fig3_cfu <- plot_cfu_data(fig3_cfudata) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(-1, 11)) 
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
#Statistical annotation labels based on adjusted kruskal-wallis p-values for first 10 days of experiment:
x_annotation <- weight_kruskal_wallis_adjust %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day)
y_position <- max(fig3_weightdata$weight_change)
label <- weight_kruskal_wallis_adjust %>% 
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>% 
  pull(p.signif)

#Narrow data to just the timepoints tested in statistical analysis, when we have weight data for all 3 groups
fig3_weightdata_subset <- fig3_weightdata %>% 
  filter(day %in% c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

fig3_weight <- plot_weight(fig3_weightdata_subset) +
  scale_x_continuous(breaks = c(-2, -4, -2, 0, 2, 4, 6, 8, 10),
                     limits = c(-3, 11)) 
save_plot(filename = "results/figures/fig3_weight.png", fig3_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

y_position <- 2 #Change for the plot showing just the median lines
fig3v2_weight <- plot_weight_medians(fig3_weightdata_subset)+
  scale_x_continuous(breaks = c(-2, -4, -2, 0, 2, 4, 6, 8, 10),
                     limits = c(-3, 11)) 
save_plot(filename = "results/figures/fig3v2_weight.png", fig3v2_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

