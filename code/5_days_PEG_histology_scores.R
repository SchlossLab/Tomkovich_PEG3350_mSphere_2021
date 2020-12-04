source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Statistical analysis and plots of histology scores from Motility #4 and #8 experiments

#Motility #4: Histology collected 6 days post-infection
d6_histology <- readxl::read_excel("data/raw/7.17.19_histo_scored_by_ingrid.xlsx", sheet = "Scored") %>%
  mutate(Tissue=factor(Tissue, levels=c("cecum", "colon")), #make sure Tissue is treated as a factor #Ignore extra column, average score for each group
         Tissue = recode(Tissue, "cecum" = "Cecum", #Capitalize cecum and colon
                                  "colon" = "Colon")) %>%
  rename(summary_score = `summary score2`) %>% #Rename summary score variable
  filter(!is.na(Group)) %>% #Drop last row, which notes that the max summary score is 12
  mutate(Group = factor(Group, levels=c("WMN", "C", "WM", "WMC")), #Change ordering of groups so that WMN (not infected) mice will appear first
         Infected = case_when(grepl("N", `Group`) ~ "no", #Make a new column based on whether mice were challenged with C. difficile
                              TRUE ~ "yes")) #An N in the Group name indicates they were not challenged with C. difficile, all the other groups were challenged
unique(d6_histology$Group) #WMN, WM, WMC, C

#Motility #8: Histology collected 4 days post-infection for most of the mice. There are 4 WMN mice that were collected at day 0
d4_histology <- readxl::read_excel("data/raw/4.5.20.histology_scored_by_ingrid.xlsx", sheet = "Scored") %>%
  mutate(Tissue=factor(Tissue, levels=c("cecum", "colon")), #make sure Tissue is treated as a factor #Ignore extra column, average score for each group
         Tissue = recode (Tissue, "cecum" = "Cecum", #Capitalize cecum and colon
                          "colon" = "Colon"),
         day = case_when(grepl("102119", `Cassette ID`) ~ "0", #Make a new column based on what date is included in casette ID
                         grepl("102519", `Cassette ID`) ~ "4"), #Tissues were collected on either day 0 or day 4
         Group= recode(Group, "WM -> WMN" = "WMN", #Recode these group names to reflect their final group identity at time of sacrifice
                              "WMN -> WM" = "WM"),  #Recode these group names to reflect their final group identity at time of sacrifice
         Group = case_when(day == 0 ~ "D0 WMN", #Rename Group to indicate these WMN mice were collected on D0 of the experiment (never challenged with C. difficile, on the last day of PEG 3350 treatment)
                           TRUE ~ Group)) %>% #Keep all other Group names the same
  rename(summary_score = `summary score (range 0-12)`) %>% #Shorten summary score column name
  mutate(Group = factor(Group, levels=c("CN", "D0 WMN", "WMN", "C", "WM")), #Change ordering of groups so that WMN (not infected) mice will appear first
         Infected = case_when(grepl("N", `Group`) ~ "no", #Make a new column based on whether mice were challenged with C. difficile
                              TRUE ~ "yes")) #An N in the Group name indicates they were not challenged with C. difficile, all the other groups were challenged

unique(d4_histology$Group) # DO WMN, WMN, CN, WM, C

#Statistical analysis: test for summary score differences across groups----
set.seed(19760620) #Same seed used for mothur analysis

#Day 6 histology analysis
#Kruskal_wallis test for differences arcoss groups in different tissue compartments
d6_kruskal_wallis_histo <- d6_histology %>%
  select(Group, summary_score, Tissue) %>%
  group_by(Tissue) %>%
  nest() %>%
  mutate(model=map(data, ~kruskal.test(x=.x$summary_score, g=as.factor(.x$Group)) %>% tidy())) %>%
  mutate(median = map(data, get_summary_score_median)) %>%
  unnest(c(model, median)) %>%
  ungroup()

#Adjust p-values for testing both cecum and colon tissues
d6_kruskal_wallis_histo_adjust <- d6_kruskal_wallis_histo %>%
  select(Tissue, statistic, p.value, parameter, method, WMN, C, WM, WMC) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj)
#Significant for colon (p = .0229), but not cecum (p = 0.105)

#Perform pairwise Wilcoxan rank sum tests for Tissues where histology scores were significantly different across groups----
d6_colon_pairwise <- d6_kruskal_wallis_histo %>%
  filter(Tissue == "Colon") %>%
  mutate(model = map(data, ~pairwise.wilcox.test(x=.x$summary_score, g=as.factor(.x$Group), p.adjust.method = "BH") %>%
                       tidy() %>%
                       mutate(compare=paste(group1, group2, sep="-")) %>%
                       select(-group1, -group2) %>%
                       pivot_wider(names_from=compare, values_from=p.value)
  )
  ) %>%
  unnest(model) %>%
  select(-data, -parameter, -statistic)

#Format pairwise stats to use with ggpubr package
d6_plot_pairwise_stats <- d6_colon_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method, -WMN, -C, -WM, -WMC) %>%
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_histology_pairwise) %>%
  bind_rows()

#Day 4 histology analysis:
d4_kruskal_wallis_histo <- d4_histology %>%
  select(Group, summary_score, Tissue) %>%
  group_by(Tissue) %>%
  nest() %>%
  mutate(model=map(data, ~kruskal.test(x=.x$summary_score, g=as.factor(.x$Group)) %>% tidy())) %>%
  mutate(median = map(data, get_summary_score_median)) %>%
  unnest(c(model, median)) %>%
  ungroup()

#Adjust p-values for testing both cecum and colon tissues
d4_kruskal_wallis_histo_adjust <- d4_kruskal_wallis_histo %>%
  select(Tissue, statistic, p.value, parameter, method, `D0 WMN`, WMN, CN, WM, C) %>%
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
  arrange(p.value.adj)
#Significant for cecum (p = .00430), but not colon (p = 0.586)

#Perform pairwise Wilcoxan rank sum tests for Tissues where histology scores were significantly different across groups----
d4_colon_pairwise <- d4_kruskal_wallis_histo %>%
  filter(Tissue == "Cecum") %>%
  mutate(model = map(data, ~pairwise.wilcox.test(x=.x$summary_score, g=as.factor(.x$Group), p.adjust.method = "BH") %>%
                       tidy() %>%
                       mutate(compare=paste(group1, group2, sep="-")) %>%
                       select(-group1, -group2) %>%
                       pivot_wider(names_from=compare, values_from=p.value)
  )
  ) %>%
  unnest(model) %>%
  select(-data, -parameter, -statistic)

#Format pairwise stats to use with ggpubr package
d4_plot_pairwise_stats <- d4_colon_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method, -`D0 WMN`, -WMN, -CN, -WM, -C) %>%
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_histology_pairwise) %>%
  bind_rows()

#Plots of summary histology scores for all mice----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMR")
color_labels <- c( "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 + 10-day recovery")
#Define shape scheme for d4 and d6 histology based on whether mice were challenged with C. difficile:
#Define shape scheme based on Infected status----
shape_scheme <- c(1, 19)
shape_infected <- c("no", "yes")

#Function to plot histology scores that were collected on day 0, 4, and 6----
#Arguments: df = dataframe to plot, pairwise_stats = dataframe of statistics to add to the plot
plot_histology <- function(df, pairwise_stats){
  group_medians <- df %>%
    group_by(Group, Tissue) %>%
    mutate(median_score = median(summary_score, na.rm = TRUE))
  histology_plot <- group_medians %>%
  ggplot(aes(x = Group, y = `summary_score`, color = Group))+
  geom_errorbar(aes(ymax = median_score, ymin = median_score), color = "gray70", size = 1.5)+ #Add lines to indicate the median for each group to the plot. Median calculated before y axis transformation
  geom_jitter(aes(shape=Infected), alpha = 0.6, size = 2)+
  labs(y = "Summary Score", x = NULL)+
  scale_color_manual(name=NULL,
                     breaks=color_groups,
                     labels=color_labels,
                     values=color_scheme)+
  scale_shape_manual(name="Infected",
                     values=shape_scheme,
                     breaks=shape_infected,
                     labels=shape_infected) +
  facet_wrap(~Tissue)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+  #So X axis labels don't overlap
  scale_y_continuous(limits = c(-0.5, 12), breaks = c(0, 4, 8, 12))+
  stat_pvalue_manual(data = pairwise_stats, label = "p.adj", y.position = "y.position", size = 6, bracket.size = .6) +
  theme_classic()+
  theme(strip.background = element_blank()) #get rid of box around facet_wrap labels
  
}

#Define color scheme for d6_histology----
color_scheme <- c("#88419d", "#238b45", "#88419d", "#f768a1") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("WMN", "C", "WM", "WMC")
color_labels <- c("5-day PEG 3350 without infection", "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.")

#Format day 6 pairwise stats to add to plot----
pairwise_day6_plot <-d6_plot_pairwise_stats %>%
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj="*") %>% #Just indicate whether statistically significant, exact p.adj values are in supplemental table
  mutate(y.position = (1:n())*4.3)

#Plot of day 6 histology scores----
d6_plot <- plot_histology(d6_histology, pairwise_day6_plot)
save_plot("results/figures/histo_scores_d6.png", d6_plot) #Use save_plot instead of ggsave because it works better with cowplot


#Plot of day 4 (plus day 0 for 1 group) histology scores----
#Define color scheme for d4_histology----
color_scheme <- c("#88419d", "#238b45", "#88419d", "#238b45", "#88419d") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("D0 WMN", "CN", "WMN", "C", "WM")
color_labels <- c("After 5-day PEG 3350 (Day 0)", "Clind. without infection", "5-day PEG 3350 without infection", "Clind.", "5-day PEG 3350")

#Format day 4 pairwise stats to add to plot----
pairwise_day4_plot <- d4_plot_pairwise_stats %>%
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj="*") %>% #Just indicate whether statistically significant, exact p.adj values are in supplemental table
  mutate(y.position = c(6, 5, 4))

#Plot of day 4 histology scores----
d4_plot <- plot_histology(d4_histology, pairwise_day4_plot)
ggsave("results/figures/histo_scores_d4.png", d4_plot)
