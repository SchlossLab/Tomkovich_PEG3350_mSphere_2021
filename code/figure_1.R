source("code/functions.R") #Reads in metadata

#Define color scheme----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMR")
color_labels <- c( "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 + 10-day recovery")

#Narrow metadata to relevant groups and experiments (WM, WMC, WMR, C)----
fig1_metadata <- metadata %>% 
  filter(group == "C" & exp_num %in% c("M3","M4", "M5", "M8")| #Only use C mice from these experiments. Allocated groups to figures based on paper outline.
         group == "WM" & exp_num %in% c("M3","M4", "M5", "M8")|
         group == "WMC" & exp_num %in% c("M3","M4")|
         group == "WMR" & exp_num %in% c("M5","M6")) %>% 
  mutate(group=factor(group, levels=c("C", "WM", "WMC", "WMR")))  # Make sure group is treated as a factor

# of mice represented in the figure
fig1_mice <- length(unique(fig1_metadata$m_id_unique)) 
# 62 mice total for figure 1

#C. difficile CFU dataframe----
#Narrow fig1_metadata to just timepoints relevant to C. difficile CFU tracking (Anything on or after day 0)
fig1_cfudata <- fig1_metadata %>% 
  filter(day > -1)
fig1_cfu_na <- sum(is.na(fig1_cfudata$avg_cfu)) #182 samples with NA values. Represent times when we either did not collect stool samples, weren't able to get a stool sample from a particular mouse, weren't able to plate the sample we did collect immediately after due to chamber issues or time constraints, or the mouse died early
#Drop rows with NA values for fig1_cfu:
fig1_cfudata <- fig1_cfudata %>% 
  filter(!is.na(avg_cfu))

#Weight change dataframe
#Note baseline weight for each group of mice (based on the earliest timepoint recorded for each experiment)----
baseline <- fig1_metadata %>% #Baseline weight was taken at day -5 for groups C, WM, and WMC
  filter(group == "C" & day == -5| #20 mice in C group
         group == "WM" & day == -5| #21 mice in WM group
         group == "WMC" & day == -5| #9 mice in WMC group
         group == "WMR" & day == -15) %>% #12 mice in WMR group, baseline weight was taken at day -15
  mutate(baseline_weight = weight) %>% #This column represents the initial weight that was recorded for each mouse
  select(m_id_unique, baseline_weight) #Will use m_id_unique to join baseline_weights to fig1_metadata

#Make a new column that represents weight_change from baseline_weight----
fig1_weightdata <- inner_join(fig1_metadata, baseline, by = "m_id_unique") %>% #Join baseline weight to fig1_metadata
  group_by(m_id_unique, day) %>% #Group by each unique mouse and experiment day
  mutate(weight_change = weight-baseline_weight) %>% #Make a new column that represents the change in weight from baseline (all weights recorded in grams)
  ungroup() %>% 
  filter(!is.na(weight)) #drop rows with NA values for fig1_weightdata. 1040 samples including NAs, 870 samples after excluding NAs

#Statistical Analysis----

#Shapiro-Wilk test to see if cfu and weight change data is normally distributed:
#Note: p-value > 0.05 means the data is normally distributed
shapiro.test(fig1_cfudata$avg_cfu) #p-value < 2.2e-16
shapiro.test(fig1_weightdata$weight_change) #p-value = 1.485e-09
#Since p-value < 0.05 for both variables, we will use non-parametric tests

#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
kruskal_wallis_cfu <- fig1_cfudata %>% 
  filter(day %in% c(0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 15, 20, 25, 30)) %>%  #only test days that we have CFU data for #Only have cfu for WMR group on D7, exclude that day
  group_by(day) %>% 
  do(tidy(kruskal.test(avg_cfu~factor(group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) 
#Timepoints where C. diff CFU is significantly different across the sources of mice
sig_C.diff_CFU_timepoints <- kruskal_wallis_cfu %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day)

#Plot data----

#Plots of CFU data----
plot_cfu <- function(df){
  median_summary <- df %>% 
    group_by(group, day) %>% 
    summarize(median_avg_cfu = median(avg_cfu, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = avg_cfu, color= group, fill = group), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(median_summary, mapping = aes(x = day, y = median_avg_cfu, color = group), alpha = 0.6, size = 1) +
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
fig1_cfu <- plot_cfu(fig1_cfudata)
save_plot(filename = "results/figures/fig1_cfu.png", fig1_cfu, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)

#plot cfu for just the inital 10days
plot_cfu_10d <- function(df){
  median_summary <- df %>% 
    group_by(group, day) %>% 
    summarize(median_avg_cfu = median(avg_cfu, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = avg_cfu, color= group, fill = group), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(median_summary, mapping = aes(x = day, y = median_avg_cfu, color = group), alpha = 0.6, size = 1) +
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
fig1_cfu_10d <- plot_cfu_10d(fig1_cfudata)
save_plot(filename = "results/figures/fig1_cfu_10d.png", fig1_cfu_10d, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)


#Weight change plot----
#Function to plot weight. Argument = dataframe you want to plot
plot_weight <- function(df){
  median_summary <- df %>% 
    group_by(group, day) %>% 
    summarize(median_weight_change = median(weight_change, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = weight_change, color= group, fill = group), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(median_summary, mapping = aes(x = day, y = median_weight_change, color = group), alpha = 0.6, size = 1) +
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

#Simplified function to plot weight that only plots the median of each group and no points for individual mice. Argument = dataframe you want to plot.
plot_weight_simple <- function(df){
  median_summary <- df %>% 
    group_by(group, day) %>% 
    summarize(median_weight_change = median(weight_change, na.rm = TRUE))
  ggplot(NULL) +
    geom_line(median_summary, mapping = aes(x = day, y = median_weight_change, color = group), alpha = 0.6, size = 1) +
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

fig1_weight <- plot_weight(fig1_weightdata) 
save_plot(filename = "results/figures/fig1_weight.png", fig1_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)
fig1v2_weight <- plot_weight_simple(fig1_weightdata)
save_plot(filename = "results/figures/fig1v2_weight.png", fig1v2_weight, base_height = 4, base_width = 8.5, base_aspect_ratio = 2)




