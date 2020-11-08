source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Reads in mothur output files

#Define color scheme to match 1 Day Peg Plots
color_scheme <- c("#225ea8", "#238b45", "#88419d") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c( "1RM1", "C", "M1")
color_labels <- c( "1-day PEG 3350 + 1-day recovery", "Clind.", "1-day PEG 3350")

# Subset alpha diversity data (16S_common_files) to analyze one day PEG subset mice
diversity_data_subset <- one_day_PEG_subset(diversity_data)
div_data_pretreatment <- diversity_data_subset %>% filter(day == -2 | day == -1)  

#Statistical Analysis----
set.seed(19760620) #Same seed used for mothur analysis

#Alpha Diversity Analysis - Shannon Plot 
diversity_data_subset <- diversity_data_subset %>%
  filter(day %in% c(-2, -1, 0, 1, 2, 4, 7))

#Plot Shannon Diversity overtime
Shannon_1_Day_PEG_Overtime <- plot_shannon_overtime(diversity_data_subset) +
  scale_x_continuous(breaks = c(-2:7),
                     limits = c(-3,8),
                     minor_breaks = c(-2.5:7.5)) +
  #Add rectangle to signify pre-treatment 
  geom_rect(mapping = aes(xmin = -2.25, xmax = -0.6, ymin = 2.5, ymax = Inf), fill = 'grey94', alpha = .01, color = 'red', show.legend = FALSE)
  save_plot("results/figures/1_Day_PEG_Overtime_shannon.png", Shannon_1_Day_PEG_Overtime)
  

#Plot sobs for 1_day_PEG subset
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
  scale_x_continuous(breaks = c(-2:7),
                     limits = c(-3, 8),
                     minor_breaks = c(-2.5:7.5))+
  theme_classic()+
  theme(legend.position = c(1,.25),
        text = element_text(size = 14), # Change font size for entire plot
        axis.ticks.x = element_blank()) +
  #Add rectange to signigy the pre-treatment
  geom_rect(mapping = aes(xmin = -2.25, xmax = -0.6, ymin = 40, ymax = Inf), fill = 'grey94', alpha = .01, color = 'red', show.legend = FALSE)
save_plot("results/figures/1_Day_PEG_richness.png", sobs_1_Day_PEG)

 #Pull 1_Day_PEG subset of PCoA data
pcoa_1_day_PEG <- read_tsv("data/process/1_day_PEG/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename("unique_label" = group) %>%
  right_join(diversity_data_subset, by= "unique_label") %>% #merge metadata and PCoA data frames (This drops some of our 16S data for early timepoints)
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


#Get exp days sequenced in both pretreatment and other days
div_data_posttreatment <- diversity_data_subset %>% filter(day >= 0)

exp_days_pretreat <- unique(div_data_pretreatment %>% pull(day))
exp_days_seq <- unique(div_data_posttreatment %>% pull(day))


#Function to test at the otu level:
kruskal_wallis_otu <- function(timepoint){
  otu_stats <- agg_otu_data %>% 
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
    write_tsv(path = paste0("data/process/otu_stats_day_", timepoint, ".tsv"))
}

# Perform kruskal wallis tests at the otu level for all days of the experiment that were sequenced----
for (d in exp_days_seq){
  kruskal_wallis_otu(d)
  #Make a list of significant otus across sources of mice for a specific day  
  stats <- read_tsv(file = paste0("data/process/otu_stats_day_", d, ".tsv"))
  name <- paste("sig_otu_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, otu))
}
  