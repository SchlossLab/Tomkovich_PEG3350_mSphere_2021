source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match 1 Day Peg Plots
color_scheme <- c("#225ea8", "#238b45", "#88419d") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c( "1RM1", "C", "M1")
color_labels <- c( "1-day PEG 3350 + 1-day recovery", "Clind.", "1-day PEG 3350")

#Subset alpha diversity data (16S_common_files) to analyze one day PEG subset mice----
diversity_data_subset <- one_day_PEG_subset(diversity_data) %>%
  mutate(day = replace(day, day == -1, "PT"),
         day = replace(day, day == -2, "PT")) #Replace day -2 and day -1 with PT to represent the pretreatment timepoint

 diversity_data_subset <- diversity_data_subset %>%
   filter(day %in% c(0, 1, 2, 4, 5, 7, "PT"))

#Get exp days sequenced in both pretreatment and other days
exp_days_seq <- unique(diversity_data_subset %>% pull(day))

#Remove Days 0 and 4 because not all three groups are represented
exp_days_seq <- exp_days_seq[1:4]  #Removes last two positions of day 0 and day 4

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Function to test at the otu level:
kruskal_wallis_otu <- function(timepoint){
  otu_stats <- agg_otu_data %>%
    mutate(day = replace(day, day == -1, "PT"),
           day = replace(day, day == -2, "PT")) %>% #Replace day -2 and day -1 with PT to represent the pretreatment timepoint
    filter(day == timepoint) %>%
    select(group, otu, agg_rel_abund) %>%
    group_by(otu) %>%
    nest() %>%
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$group)) %>% tidy())) %>%
    mutate(median = map(data, get_rel_abund_median_group)) %>%
    unnest(c(model, median)) %>%
    ungroup()
  #Adjust p-values for testing multiple OTUs
  otu_stats_adjust <- otu_stats %>%
    select(otu, statistic, p.value, parameter, method, "C", "1RM1", "M1") %>%
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>%
    arrange(p.value.adj) %>%
    write_tsv(path = paste0("data/process/1_Day_PEG_otu_stats_day_", timepoint, ".tsv"))
}

# Perform kruskal wallis tests at the otu level for all days of the experiment that were sequenced----
for (d in exp_days_seq){
  kruskal_wallis_otu(d)
  #Make a list of significant otus across sources of mice for a specific day
  stats <- read_tsv(file = paste0("data/process/1_Day_PEG_otu_stats_day_", d, ".tsv"))
  name <- paste("sig_otu_day", d, sep = "")
  assign(name, pull_significant_taxa(stats, otu))
}

#Shared significant genera across from Days 1, 2 and 7----
shared_sig_otus_D1toD7 <- intersect_all(`sig_otu_day1`, sig_otu_day2, sig_otu_day5, sig_otu_day7, sig_otu_dayPT)
view(shared_sig_otus_D1toD7)
print(shared_sig_otus_D1toD7)
#[31] "Bacteroides (OTU 1)"         "Enterobacteriaceae (OTU 2)"  "Porphyromonadaceae (OTU 8)" 
#[4] "Lachnospiraceae (OTU 20)"    "Lachnospiraceae (OTU 62)"    "Turicibacter (OTU 5)"       
#[7] "Porphyromonadaceae (OTU 14)" "Porphyromonadaceae (OTU 6)"  "Porphyromonadaceae (OTU 7)" 
#[10] "Lachnospiraceae (OTU 32)"    "Lachnospiraceae (OTU 4)"     "Porphyromonadaceae (OTU 10)"
#[13] "Lachnospiraceae (OTU 11)"    "Lachnospiraceae (OTU 24)"    "Akkermansia (OTU 3)"    

#Function to plot a list of OTUs across sources of mice at a specific timepoint:
#Arguments: otus = list of otus to plot; timepoint = day of the experiment to plot
plot_otus_dx <- function(otus, timepoint){
  one_day_PEG_subset(agg_otu_data) %>%
    mutate(day = replace(day, day == -1, "PT"),
           day = replace(day, day == -2, "PT")) %>% #Replace day -2 and day -1 with PT to represent the pretreatment timepoint
    filter(otu %in% otus) %>%
    filter(day == timepoint) %>%
    mutate(agg_rel_abund = agg_rel_abund + 1/2000) %>% # 2000 is 2 times the subsampling parameter of 1000
    ggplot(aes(x= otu_name, y=agg_rel_abund, color=group))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
    geom_hline(yintercept=1/2000, color="gray")+
    stat_summary(fun = 'median',
                 fun.max = function(x) quantile(x, 0.75),
                 fun.min = function(x) quantile(x, 0.25),
                 position = position_dodge(width = 1)) +
    labs(title=NULL,
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/2000, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "bottom",
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Plots of the top 20 OTUs that varied across sources at each timepoint----
#Day 1: top 20 OTUs that vary across sources
D1top20_otus <- plot_otus_dx(`sig_otu_day1`[1:20], 1) +#Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_D1top20_otus.png", D1top20_otus, base_height = 9, base_width = 7)
#Day 2: top 20 OTUs that vary across sources
D2top20_otus <- plot_otus_dx(`sig_otu_day2`[1:20], 2) + #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_D2top20_otus.png", D2top20_otus, base_height = 9, base_width = 7)
#Day 5: top 20 OTUs that vary across sources
D5top20_otus <- plot_otus_dx(`sig_otu_day5`[1:20], 5)+ #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_D5top20_otus.png", D5top20_otus, base_height = 9, base_width = 7)
#Day 7: top 20 OTUs that vary across sources
D7top20_otus <- plot_otus_dx(`sig_otu_day7`[1:20], 7)+ #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_D7top20_otus.png", D7top20_otus, base_height = 9, base_width = 7)
#Pre-treatment (PT): top 20 OTUs that vary across sources
PTtop20_otus <- plot_otus_dx(`sig_otu_dayPT`[1:20], "PT")+ #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/1_Day_PEG_PTtop20_otus.png", PTtop20_otus, base_height = 9, base_width = 7)

#Plot sobs overtime for 1_day_PEG subset----
sobs_1_Day_PEG <- diversity_data_subset %>%
  group_by(group, day) %>%
  mutate(median_sobs = median(sobs)) %>%
  ggplot(x = day, y = sobs, colour =  group) +
  geom_point(mapping = aes(x = day, y = sobs, color = group, fill = group), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
  geom_line(mapping = aes(x = day, y = median_sobs, color = group), alpha = 0.6, size = 1) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels) +
  theme_classic()+
  theme(legend.position = c(1,.25),
        text = element_text(size = 14), # Change font size for entire plot
        axis.ticks.x = element_blank()) 

save_plot("results/figures/1_Day_PEG_sobs_overtime.png", sobs_1_Day_PEG)

#Plot Shannon Diversity overtime
Shannon_1_Day_PEG_Overtime <- plot_shannon_overtime(diversity_data_subset) +
  scale_x_continuous(breaks = c(-2:7),
                     limits = c(-3,8),
                     minor_breaks = c(-2.5:7.5)) +
  #Add rectangle to signify pre-treatment
  geom_rect(mapping = aes(xmin = -2.25, xmax = -0.6, ymin = 2.5, ymax = Inf), fill = 'grey94', alpha = .01, color = 'red', show.legend = FALSE)
save_plot("results/figures/1_Day_PEG__shannon_overtime.png", Shannon_1_Day_PEG_Overtime)

 #Pull 1_Day_PEG subset of PCoA data
pcoa_1_day_PEG <- read_tsv("data/process/1_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  full_join(diversity_data_subset, by= "unique_label") %>% #merge metadata and PCoA data frames (This drops some of our 16S data for early timepoints)
  mutate(day = as.integer(day)) %>% #Day variable (transformed to integer to get rid of decimals on PCoA animation
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff


#Pull axes from loadings file
pcoa_axes_1_day_PEG <- read_tsv("data/process/1_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
axis1 <- pcoa_axes_1_day_PEG %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
axis2 <- pcoa_axes_1_day_PEG %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal


pcoa_subset_plot <- plot_pcoa(pcoa_1_day_PEG)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))
save_plot(filename = paste0("results/figures/1_Day_PEG_PCoA.png"), pcoa_subset_plot, base_height = 5, base_width = 4.5)

pcoa_animated <- plot_pcoa(pcoa_1_day_PEG)+
  labs(x = paste("PCoA 1 (", axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", axis2,"%)", sep = ""))+
  labs(title = 'Day: {frame_time}') + #Adds time variable to title
  transition_time(day)+  #Day variable used to cycle through time on animation
  shadow_mark() #Shows previous timepoints

# Implement better frames per second for animation
pcoa_gif <- animate(pcoa_animated, duration = 7, fps = 10,
                    res = 150, width = 20, height = 20, unit = "cm")

# Save as gif file
anim_save(animation = pcoa_gif, filename = 'results/1_Day_PEG_pcoa_over_time.gif')
