source("code/functions.R") #Loads libraries, reads in metadata, functions

#Define color scheme to match figure 1 plots----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#225ea8") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMR")
color_labels <- c( "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 + 10-day recovery")

metadata <- metadata %>% 
  mutate(day = as.integer(day)) %>%  #Day variable (transformed to integer to get rid of decimals on PCoA animation
  select(-stool_tube_label, -tissue_type, -tissue_tube_label) %>% #Get rid of columns not needed
  mutate(day = as.integer(day))

#Lost 3 samples since we subsampled to 2000 sequences per sample
pcoa_data <- read_tsv("data/process/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(unique_label = group) %>% #group is the same as id in the metadata data frame
  mutate(unique_label = replace(unique_label, unique_label == "M5WMR10D11", "M5WMR10D1"), #correct unique_labels that have typos so that they'll match up with metadata
         unique_label = replace(unique_label, unique_label == "M5WMR5D11", "M5WMR5D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR5D14", "M5WMR5D4"), 
         unique_label = replace(unique_label, unique_label == "M5WMR6D11", "M5WMR6D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR6D14", "M5WMR6D4"),
         unique_label = replace(unique_label, unique_label == "M5WMR7D11", "M5WMR7D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR7D14", "M5WMR7D4"),
         unique_label = replace(unique_label, unique_label == "M5WMR8D11", "M5WMR8D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR8D14", "M5WMR8D4"),
         unique_label = replace(unique_label, unique_label == "M5WMR9D11", "M5WMR9D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR9D14", "M5WMR9D4")) %>% 
  filter(unique_label != "M8C18D4") %>% #Remove this sample, since I discovered this mouse was pregnant during tissue collection
  left_join(metadata, by= "unique_label") #merge metadata and PCoA data frames

test <- pcoa_data %>% count(day)

#Statistical Analysis----
set.seed(19881117) #Match seed used in mothur analysis scripts

fig1_dist <- read_dist("data/process/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")

#Get factor levels for mouse_id variable:
mouse_id_levels <- unique(as.factor(pcoa_data$m_id_unique))
#48 levels

#Get factor levels for unique_cage variable:
unique_cage_levels <- unique(as.factor(pcoa_data$unique_cage))
#25 levels

#Plot PCoA data----
#Function to plot PCoA data 
plot_pcoa <- function(df){
  ggplot(df, aes(x=axis1, y=axis2, color = group, alpha = day)) +
    geom_point(size=2) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_groups,
                        labels=color_labels)+
#    scale_alpha_continuous(range = c(.3, 1),
#                           breaks= c(2, 4, 6, 8, 10),
#                           labels=c(2, 4, 6, 8, 10))+
    coord_fixed() + 
#    xlim(-0.4, 0.65)+
#    ylim(-0.45, 0.6)+
    labs(x="PCoA 1",
         y="PCoA 2",
         alpha= "Day") +
    theme_classic()
}

#Read in .loadings file to add percent variation represented by PCoA axis
all_axis_labels <- read_tsv("data/process/peg3350.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings")
all_axis1 <- all_axis_labels %>% filter(axis == 1) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal
all_axis2 <- all_axis_labels %>% filter(axis == 2) %>% pull(loading) %>% round(digits = 1) #Pull value & round to 1 decimal

#PCoA plot that combines the 2 experiments and save the plot----  
fig1_pcoa_plot <- plot_pcoa(pcoa_data)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))
save_plot(filename = paste0("results/figures/fig1_pcoa.png"), fig1_pcoa_plot, base_height = 5, base_width = 4.5)

fig1_pcoa_plot_time <- plot_pcoa(pcoa_data)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))+
  theme(legend.position = "none")+ #remove legend
  facet_wrap(~ day)

fig1_pcoa_plot_dn5 <- plot_pcoa(pcoa_data %>% filter(day == "-5"))+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))

fig1_pcoa_plot_d1 <- plot_pcoa(pcoa_data %>% filter(day == "1"))+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))

fig1_pcoa_plot_d4 <- plot_pcoa(pcoa_data %>% filter(day == "4"))+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))

fig1_pcoa_plot_d10 <- plot_pcoa(pcoa_data %>% filter(day == "10"))+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))

#Animation of PCoA plot over time for all sequenced samples ----
#Source: Will Close's Code Club from 4/12/2020 on plot animation
pcoa_animated <- plot_pcoa(pcoa_data)+
  labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Anotations for each axis from loadings file
       y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))+
  labs(title = 'Day: {frame_time}') + #Adds time variable to title
  transition_time(day)+  #Day variable used to cycle through time on animation
  shadow_mark() #Shows previous timepoints

# Implement better frames per second for animation
pcoa_gif <- animate(pcoa_animated, duration = 10, fps = 10,
                    res = 150, width = 20, height = 20, unit = "cm")

# Save as gif file
anim_save(animation = pcoa_gif, filename = 'results/fig1_pcoa_over_time.gif')

#Alpha diversity analysis----
diversity_data <- read_tsv("data/process/peg3350.opti_mcc.groups.ave-std.summary") %>%
  filter(method == "ave") %>% 
  select(group, sobs, shannon, invsimpson, coverage) %>%   
  rename(unique_label = group) %>% #group is the same as unique_label in the metadata data frame
  mutate(unique_label = replace(unique_label, unique_label == "M5WMR10D11", "M5WMR10D1"), #correct unique_labels that have typos so that they'll match up with metadata
         unique_label = replace(unique_label, unique_label == "M5WMR5D11", "M5WMR5D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR5D14", "M5WMR5D4"), 
         unique_label = replace(unique_label, unique_label == "M5WMR6D11", "M5WMR6D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR6D14", "M5WMR6D4"),
         unique_label = replace(unique_label, unique_label == "M5WMR7D11", "M5WMR7D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR7D14", "M5WMR7D4"),
         unique_label = replace(unique_label, unique_label == "M5WMR8D11", "M5WMR8D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR8D14", "M5WMR8D4"),
         unique_label = replace(unique_label, unique_label == "M5WMR9D11", "M5WMR9D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR9D14", "M5WMR9D4")) %>% 
  filter(unique_label != "M8C18D4") %>% #Remove this sample, since I discovered this mouse was pregnant during tissue collection
  inner_join(metadata, by = "unique_label") #Match only the samples we have sequence data for

#Plot of shannon diversity at days 1, 4, and 10 when we have sequencing data for 3 groups
shannon_d1_4_10 <- diversity_data %>% 
  filter(day %in% c(1, 4, 10)) %>% 
  group_by(group, day) %>% 
  mutate(median_shannon = median(shannon)) %>% #create a column of median values for each group
  ungroup() %>% 
  ggplot(aes(x=group, y =shannon, colour= group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
#  scale_shape_manual(name=NULL,
#                     values=shape_scheme,
#                     breaks=shape_experiment,
#                     labels=shape_experiment) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median_shannon, ymin = median_shannon), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(size=2, alpha=0.6, show.legend = FALSE) +
  labs(title=NULL, 
       x=NULL,
       y="Shannon Diversity Index")+
  ylim(0, 3.5)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank()) #Remove x axis ticks
save_plot("results/figures/shannon_d1_4_10.png", shannon_d1_4_10) #Use save_plot instead of ggsave because it works better with cowplot

#Plot of WMR group over the days we have sequencing data for 3 groups
shannon_WMR <- diversity_data %>% 
  filter(group == "WMR") %>% 
  group_by(day) %>% 
  mutate(median_shannon = median(shannon)) %>% #create a column of median values for each group
  ungroup() %>% 
  ggplot(aes(x=group, y =shannon, colour= group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  #  scale_shape_manual(name=NULL,
  #                     values=shape_scheme,
  #                     breaks=shape_experiment,
  #                     labels=shape_experiment) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median_shannon, ymin = median_shannon), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(size=2, alpha=0.6, show.legend = FALSE) +
  labs(title=NULL, 
       x=NULL,
       y="Shannon Diversity Index")+
  ylim(0, 3.5)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank()) #Remove x axis ticks
save_plot("results/figures/shannon_WMR.png", shannon_WMR) #Use save_plot instead of ggsave because it works better with cowplot

#Plot of sobs (richness) at days 1, 4, and 10 when we have sequencing data for 3 groups
sobs_d1_4_10 <- diversity_data %>% 
  filter(day %in% c(1, 4, 10)) %>% 
  group_by(group, day) %>% 
  mutate(median_sobs = median(sobs)) %>% #create a column of median values for each group
  ungroup() %>% 
  ggplot(aes(x=group, y =sobs, colour= group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  #  scale_shape_manual(name=NULL,
  #                     values=shape_scheme,
  #                     breaks=shape_experiment,
  #                     labels=shape_experiment) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median_sobs, ymin = median_sobs), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(size=2, alpha=0.6, show.legend = FALSE) +
  labs(title=NULL, 
       x=NULL,
       y="Number of Observed OTUs")+
  ylim(0, 150)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank()) #Remove x axis ticks
save_plot("results/figures/richness_d1_4_10.png", sobs_d1_4_10) #Use save_plot instead of ggsave because it works better with cowplot

#Plot of sobs (richness) for the WMR group over the days we have sequencing data for 3 groups
sobs_WMR <- diversity_data %>% 
  filter(group == "WMR") %>% 
  group_by(day) %>% 
  mutate(median_sobs = median(sobs)) %>% #create a column of median values for each group
  ungroup() %>% 
  ggplot(aes(x=group, y =sobs, colour= group))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_groups,
                      labels=color_labels)+
  #  scale_shape_manual(name=NULL,
  #                     values=shape_scheme,
  #                     breaks=shape_experiment,
  #                     labels=shape_experiment) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median_sobs, ymin = median_sobs), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(size=2, alpha=0.6, show.legend = FALSE) +
  labs(title=NULL, 
       x=NULL,
       y="Number of Observed OTUs")+
  ylim(0, 150)+
  facet_wrap(~ day)+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16), # Change font size for entire plot
        axis.text.x= element_blank(),#Remove x axis labels
        axis.ticks.x = element_blank()) #Remove x axis ticks
save_plot("results/figures/richness_WMR.png", sobs_WMR) #Use save_plot instead of ggsave because it works better with cowplot

#OTU analysis----
# Import taxonomy into data frame and clean up taxonomy names
taxonomy <- read_tsv(file="data/process/peg3350.taxonomy") %>%
  rename_all(tolower) %>% #remove uppercase from column names
  # Split taxonomic information into separate columns for each taxonomic level  
  mutate(taxonomy=str_replace_all(taxonomy, c("\\(\\d*\\)" = "", #drop digits with parentheses around them
                                              ';$' = "", #removes semi-colon at end of line
                                              'Bacteria_unclassified' = 'Unclassified',
                                              "Clostridium_" = "Clostridium ", #Remove underscores after Clostridium
                                              "_" = " ", #Removes all other underscores
                                              "unclassified" = "Unclassified"))) %>% 
  # Separate taxonomic levels into separate columns according to semi-colon.
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';')
#Warning message: ~80 OTUs that are unknown, filled in with NAs

# Import otu_data for samples
otu_data <- read_tsv("data/process/peg3350.opti_mcc.shared", col_types=cols(Group=col_character())) %>% 
  select(-label, -numOtus) %>% 
  rename(unique_label = Group) %>% #group is the same as unique_label in the metadata data frame
  mutate(unique_label = replace(unique_label, unique_label == "M5WMR10D11", "M5WMR10D1"), #correct unique_labels that have typos so that they'll match up with metadata
         unique_label = replace(unique_label, unique_label == "M5WMR5D11", "M5WMR5D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR5D14", "M5WMR5D4"), 
         unique_label = replace(unique_label, unique_label == "M5WMR6D11", "M5WMR6D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR6D14", "M5WMR6D4"),
         unique_label = replace(unique_label, unique_label == "M5WMR7D11", "M5WMR7D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR7D14", "M5WMR7D4"),
         unique_label = replace(unique_label, unique_label == "M5WMR8D11", "M5WMR8D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR8D14", "M5WMR8D4"),
         unique_label = replace(unique_label, unique_label == "M5WMR9D11", "M5WMR9D1"),
         unique_label = replace(unique_label, unique_label == "M5WMR9D14", "M5WMR9D4")) %>% 
  filter(unique_label != "M8C18D4") %>% #Remove this sample, since I discovered this mouse was pregnant during tissue collection
  gather(-unique_label, key="otu", value="count") %>% 
  mutate(rel_abund=count/2000) #Use 2000, because this is the subsampling parameter chosen.

#Merge otu_data to taxonomy data frame
agg_taxa_data <- inner_join(otu_data, taxonomy)

# Function to summarize relative abundance level for a given taxonomic level (ex. genus, family, phlyum, etc.)
agg_taxonomic_data <- function(taxonomic_level) {
  agg_taxa_data %>% 
    group_by(unique_label, {{ taxonomic_level }}) %>% #Embracing treats the taxonomic_level argument as a column name
    summarize(agg_rel_abund=sum(rel_abund)) %>% 
    # Merge relative abundance data to specifci taxonomic_level data
    inner_join(., metadata, by = "unique_label") %>% 
    ungroup() 
}

# Relative abundance data at the otu level:
agg_otu_data <- agg_taxonomic_data(otu)

#Rename otus to match naming convention used previously and add a column that will work with ggtext package:
agg_otu <- agg_otu_data %>% 
  mutate(key=otu) %>% 
  group_by(key)
taxa_info <- read.delim('data/process/peg3350.taxonomy', header=T, sep='\t') %>% 
  select(-Size) %>% 
  mutate(key=OTU) %>% 
  select(-OTU)
agg_otu_data <- inner_join(agg_otu, taxa_info, by="key") %>%
  ungroup() %>% 
  mutate(key=str_to_upper(key)) %>% 
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>% 
  mutate(taxa=gsub("(.*)_.*","\\1",Taxonomy)) %>% 
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>% 
  mutate(taxa=gsub(".*;","",taxa)) %>% 
  mutate(taxa=gsub("(.*)_.*","\\1",taxa)) %>% 
  mutate(taxa=gsub('[0-9]+', '', taxa)) %>% 
  mutate(taxa=str_remove_all(taxa, "[(100)]")) %>% 
  unite(key, taxa, key, sep=" (") %>% 
  mutate(key = paste(key,")", sep="")) %>% 
  select(-otu, -Taxonomy) %>% 
  rename(otu=key) %>% 
  mutate(otu=paste0(gsub('TU0*', 'TU ', otu))) %>% 
  separate(otu, into = c("bactname", "OTUnumber"), sep = "\\ [(]", remove = FALSE) %>% #Add columns to separate bacteria name from OTU number to utilize ggtext so that only bacteria name is italicized
  mutate(otu_name = glue("*{bactname}* ({OTUnumber}")) #Markdown notation so that only bacteria name is italicized

#Figure out which days we have sequencing data for from the 3 groups:
test <- pcoa_data %>% group_by(group) %>% count(day)
#Days that we have data for all 3 groups: 1, 4, 10
test_days <- c(1, 4, 10)

#Function to test for differences across groups at the OTU level for specific timepoints
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
    select(otu, statistic, p.value, parameter, method, C, WM, WMR) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/fig1_otu_stats_day_", timepoint, ".tsv"))
}

# Perform kruskal wallis tests at the otu level for all days of the experiment that were sequenced----
for (d in test_days){
  kruskal_wallis_otu(d)
  #Make a list of significant otus across sources of mice for a specific day  
  stats <- read_tsv(file = paste0("data/process/fig1_otu_stats_day_", d, ".tsv"))
  name <- paste("sig_otu_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, otu))
}

#OTUs that varied across treatment groups and were shared across days 1, 4 and 10
shared_sig_otus_d1_4_10 <- intersect_all(sig_otu_day1, sig_otu_day4, sig_otu_day10)
shared_sig_otus_d1_4 <- intersect_all(sig_otu_day1, sig_otu_day4)
shared_sig_otus_d1_10 <- intersect_all(sig_otu_day1, sig_otu_day10)
shared_sig_otus_d4_10 <- intersect_all(sig_otu_day4, sig_otu_day10)

#Only Lachnospiraceae (OTU 134) varies across all 3 timepoints
#7 OTUs vary on day 1: 1 OTU shared with day 10
#124 OTUs on day 4: 6 OTus shared with day 1
#25 OTUs on day 10: 23/25 shared with day 4

#Function to plot a list of OTUs across sources of mice at a specific timepoint:
#Arguments: otus = list of otus to plot; timepoint = day of the experiment to plot
plot_otus_dx <- function(otus, timepoint){
  agg_otu_data %>% 
    filter(otu %in% otus) %>% 
    filter(day == timepoint) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/4000) %>% # 4,000 is 2 times the subsampling parameter of 2000
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
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "none",
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Plots of the relative abundances of OTUs that significantly varied across sources of mice from day -1 to day 1----
otus_d1 <- plot_otus_dx(sig_otu_day1, 1)+
  ggtitle("Day 1 post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/fig1_otus_d1.png", otus_d1, base_height = 7, base_width = 8)

otus_d4 <- plot_otus_dx(sig_otu_day4, 4)+
  ggtitle("Day 4 post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/fig1_otus_d4.png", otus_d4, base_height = 7, base_width = 8)

otus_d10 <- plot_otus_dx(sig_otu_day10, 10)+
  ggtitle("Day 10 post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/fig1_otus_d10.png", otus_d10, base_height = 7, base_width = 8)

#Examine impacts of clindamycin and PEG3350 treatments on bacterial OTUs----

#Examine changes that happen after clindamycin treatment (baseline day -5 versus day 1)
C_dn5_d1_pairs <- agg_otu_data %>% 
  filter(group == "C" & otu == "Bacteroides (OTU 1)") %>% #Limit to group "C" and randomly pick an OTU just to figure out what mice have sequence data
  filter(day == -5 | day == 1) %>% 
  filter(duplicated(m_id_unique)) %>% #Pull mouse ids with sequence data for both day -1 and day 0
  pull(m_id_unique) #6 mice

#Dataframe for statistical test at the OTU level
C_paired_otu <- agg_otu_data %>% 
  filter(m_id_unique %in% C_dn5_d1_pairs) %>% #Only select pairs with data for day -1 & day 0
  filter(day == -5 | day == 1) %>% #Experiment days that represent initial community and community post clindamycin treatment
  mutate(day = as.factor(day)) %>% 
  select(day, otu, agg_rel_abund)

#Wilcoxon signed rank test for all day -1, day 0 pairs at the OTU level:
otus_C_pairs <- C_paired_otu %>% 
  group_by(otu) %>% 
  nest() %>% 
  mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>% 
  mutate(median = map(data, get_rel_abund_median_day)) %>% 
  unnest(c(model, median)) %>% 
  ungroup() 

#Adjust p-values for testing multiple OTUs
otus_C_pairs_stats_adjust <- otus_C_pairs %>% 
  select(otu, statistic, p.value, method, alternative, `-5`, `1`) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>% 
  write_tsv(path = "data/process/fig1_otu_Cgroup_dn5to0.tsv") 

#Make a list of significant OTUs impacted by clindamycin treatment----  
C_sig_otu_pairs <- pull_significant_taxa(otus_C_pairs_stats_adjust, otu)
# 0 OTUs
C_sig_otu_pairs_top10 <- C_sig_otu_pairs[1:10]

C_top_OTUs <- otus_C_pairs_stats_adjust %>% 
  arrange(p.value)
C_top10_OTUs <- head(C_top_OTUs, 10) %>% pull(otu)

#Examine OTUs that change after 5-day PEG3350 treatment (baseline day -5 versus day 1)
WM_dn5_d1_pairs <- agg_otu_data %>% 
  filter(group == "WM" & otu == "Bacteroides (OTU 1)") %>% #Limit to group "C" and randomly pick an OTU just to figure out what mice have sequence data
  filter(day == -5 | day == 1) %>% 
  filter(duplicated(m_id_unique)) %>% #Pull mouse ids with sequence data for both day -1 and day 0
  pull(m_id_unique) #9 mice

#Dataframe for statistical test at the OTU level
WM_paired_otu <- agg_otu_data %>% 
  filter(m_id_unique %in% WM_dn5_d1_pairs) %>% #Only select pairs with data for day -1 & day 0
  filter(day == -5 | day == 1) %>% #Experiment days that represent initial community and community post clindamycin treatment
  mutate(day = as.factor(day)) %>% 
  select(day, otu, agg_rel_abund)

#Wilcoxon signed rank test for all day -1, day 0 pairs at the OTU level:
otus_WM_pairs <- WM_paired_otu %>% 
  group_by(otu) %>% 
  nest() %>% 
  mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>% 
  mutate(median = map(data, get_rel_abund_median_day)) %>% 
  unnest(c(model, median)) %>% 
  ungroup() 

#Adjust p-values for testing multiple OTUs
otus_WM_pairs_stats_adjust <- otus_WM_pairs %>% 
  select(otu, statistic, p.value, method, alternative, `-5`, `1`) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>% 
  write_tsv(path = "data/process/fig1_otu_WMgroup_dn5to0.tsv") 

#Make a list of significant OTUs impacted by PEG3350 treatment----  
WM_sig_otu_pairs <- pull_significant_taxa(otus_WM_pairs_stats_adjust, otu)
# 0 OTUs
WM_sig_otu_pairs_top10 <- WM_sig_otu_pairs[1:10]

WM_top_OTUs <- otus_WM_pairs_stats_adjust %>% 
  arrange(p.value)
WM_top10_OTUs <- head(WM_top_OTUs, 10) %>% pull(otu)

#Compare OTUs impacted by clindamycin and PEG3350
WM_C_otus <- intersect_all(C_top10_OTUs, WM_top10_OTUs)
#2 OTUs overlap: Enterobacteriaceae (OTU 3) and Porphyromonadaceae (OTU 16)


