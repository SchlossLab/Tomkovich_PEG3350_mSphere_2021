source("code/utilities.R") #Loads libraries, reads in metadata, functions
source("code/16S_common_files.R") #Loads 16S output files for all sequenced samples

#Examine the number of sequences per sample to determine subsampling parameter----

data <- read_tsv("data/mothur/peg3350.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count.summary", col_names=c("unique_label", "nseqs")) %>% 
  left_join(metadata, by = "unique_label")
#1509 samples

#Range of sequences per sample in water controls:
waters <- data %>% 
  filter(str_detect(unique_label, "water")) %>% 
  summary(nseqs)
# Range: 4-5121 with median = 107

#Range of sequences per sample in FMT & PBS gavages:
gavage_sol <- data %>% 
  filter(str_detect(unique_label, "FMT") | str_detect(unique_label, "PBS")) %>% 
  summary(nseqs)
# Range: 197-12240 with median = 6911. PBS was the sample with 197 sequences

#Range of sequences per sample in the tissue samples:
tissues <- data %>% 
  filter(sample_type == "distal_colon"|sample_type == "proximal_colon"|sample_type == "cecum")%>% 
  summary(nseqs)
# Range: 18-94104  with median = 17508

#Range of sequences per sample for rest of stool samples:
stools <- data %>%   
  filter(!str_detect(unique_label, "water")) %>% 
  filter(!str_detect(unique_label, "FMT")) %>% 
  filter(!str_detect(unique_label, "PBS")) %>%   
  filter(!sample_type == "distal_colon") %>% 
  filter(!sample_type == "proximal_colon") %>% 
  filter(!sample_type == "cecum") %>% 
  summary(nseqs)
# Range: 9-166261  with median = 17848

#Examine number of sequences per sample----
data %>% ggplot(aes(x=nseqs))+ geom_histogram() + theme_classic() 

data %>% ggplot(aes(x=nseqs)) + geom_histogram() + scale_x_log10(limits = c(-1, 200050)) +
  geom_vline(xintercept = 1000, linetype = 2, color = "grey40")+ #Add line to note 5000 sequences per sample threshold 
  theme_classic()+
  theme(text = element_text(size = 16))+  # Change font size for entire plot
  ggsave("exploratory/notebook/seq_per_sample_distribution.pdf")

#Explore number of samples that will be lost depending on subsampling parameter chosen----

#Rarefy to 1000:
n_1000 <- data %>% filter(nseqs < 1000) %>% select(unique_label) %>% nrow()
#Lose 125 samples (number includes 14 water controls & PBS gavage solution)

#Rarefy to 500:
n_500 <- data %>% filter(nseqs < 500) %>% select(unique_label) %>% nrow()
#Lose 89 samples

#Rarefy to 1500:
n_1500 <- data %>% filter(nseqs < 1500) %>% select(unique_label) %>% nrow()
#Lose 158 samples.

#Rarefy to 2000:
n_2000 <- data %>% filter(nseqs < 2000) %>% select(unique_label) %>% nrow()
#Lose 189 samples.

#Rarefy to 4000:
n_4000 <- data %>% filter(nseqs < 4000) %>% select(unique_label) %>% nrow()
#Lose 291 samples.

# Examination of the OTUs in the water controls with > 1000 sequences----
#Figure out the top 10 OTUs across the 3 water controls:
top_10_water_otus <- agg_otu_data %>% 
  filter(str_detect(unique_label, "water")) %>% 
  mutate(sample_type = "water") %>% 
  group_by(sample_type, otu) %>% 
  summarise(mean_rel_abund = mean(agg_rel_abund)) %>% 
  top_n(10, mean_rel_abund) %>% 
  pull(otu)
#Plot the top 10 OTUs across the 3 water controls:
water_otus <- agg_otu_data %>% 
  filter(str_detect(unique_label, "water")) %>% #select only the water controls
  filter(otu %in% top_10_water_otus) %>% #Select the top 10 water OTUs
  ggplot(aes(fill = otu, y = agg_rel_abund, x = unique_label))+
  geom_bar(position = "fill", stat="identity")+
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance")+
  coord_flip()+
  theme_classic()+
  ggsave("exploratory/notebook/waters_top_10_otus.pdf")

#Figure out the top 10 OTUs in each of the 3 types of tissues: 
top_10_tissues <- agg_otu_data %>% 
  filter(!str_detect(unique_label, "water")) %>% #Remove water controls
  filter(!str_detect(unique_label, "FMT")) %>% #Remove FMT gavage samples
  filter(sample_type != "stool") %>% #Remove stool samples
  group_by(sample_type, otu) %>% 
  summarise(mean_rel_abund = mean(agg_rel_abund)) %>% 
  top_n(10, mean_rel_abund) %>% 
  pull(otu)
#Plot the top 10 OTUs in each type of tissue
sample_type_otus <- agg_otu_data %>% 
  #Replace NAs for any of the relevant variables (primarily just for the water and FMT samples)
  mutate(sample_type = case_when(str_detect(unique_label, "water") ~ "water", 
                                 str_detect(unique_label, "FMT") ~ "FMT",
                                 TRUE ~ sample_type)) %>% 
  filter(!str_detect(unique_label, "water")) %>% #Remove water controls or leave the water controls to see how they compare to the rest of the tissues
  filter(!str_detect(unique_label, "FMT")) %>% #Remove FMT gavage samples
  filter(sample_type != "stool") %>% #Remove stool samples
  filter(otu %in% top_10_tissues) %>% #Select the top 10 water OTUs
  ggplot(aes(fill = otu, y = agg_rel_abund, x = sample_type))+
  geom_bar(position = "fill", stat="identity")+
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance")+
  coord_flip()+
  theme_classic()+
  ggsave("exploratory/notebook/tissues_top_10_otus.pdf")
  
# See how plate number and library number contribute to community structure----
#Statistical analysis of PCoA data (all samples)
set.seed(19760620) #Same seed used for mothur analysis
# Read in Bray-Curtis distance matrix that represents all sequenced samples
all_dist <- read_dist("data/process/peg3350.opti_mcc.braycurtis.0.03.lt.ave.dist")

#Get factor levels for unique_mouse_id variable:
mouse_id_levels <- unique(as.factor(metadata$unique_mouse_id))
#Get factor levels for unique_cage_no variable:
cage_levels <- unique(as.factor(metadata$unique_cage_no))
#Get factor levels for group variable:
group_levels <- unique(as.factor(metadata$group))
#Get factor levels for exp_num variable:
exp_levels <- unique(as.factor(metadata$exp_num))
#Get factor levels for day variable:
day_levels <- unique(as.factor(metadata$day))
#Get factor levels for sample_type variable
sample_levels <- unique(as.factor(metadata$sample_type))

#Select relevant variables from seq_prep_metadata (the plate and MiSeq library numbers the sample was sequenced in)
prep_variables <- seq_prep_metadata %>% 
  select(unique_label, ext_plate, miseq_run) 
  
#Remove the 2nd sequencing entry for the few samples that were sequenced twice 
prep_remove <- seq_prep_metadata %>% 
  filter(unique_label == "M3WM1Dn5" & miseq_run == "CDIp52_Motility11-13"|unique_label == "M3WM3D1" & miseq_run == "CDIp52_Motility11-13"|unique_label == "M5WMR5D10" & miseq_run == "CDIp52_Motility11-13"|unique_label == "M5WMR7D10" & miseq_run == "CDIp52_Motility11-13")

#Remove the 4 samples with duplicate entries:
prep_variables <- anti_join(prep_variables, prep_remove)

#Extract sample ids from distance matrix and join to metadata and prep_variables in order to test the impact of the relevant variables on community structure
all_variables <- tibble(unique_label = attr(all_dist, "Labels")) %>% 
  left_join(metadata, by = "unique_label") %>% 
  left_join(prep_variables, by = "unique_label") %>% 
  #Replace NAs for any of the relevant variables (primarily just for the water and FMT samples)
  mutate(sample_type = case_when(str_detect(unique_label, "water") ~ "water",
                                 str_detect(unique_label, "FMT") ~ "FMT",
                                 TRUE ~ sample_type)) %>% 
  mutate(unique_mouse_id = replace_na(unique_mouse_id, "not_applicable")) %>% 
  mutate(unique_cage_no = replace_na(unique_cage_no, "not_applicable")) %>% 
  mutate(group = replace_na(group, "not_applicable")) %>% 
  mutate(exp_num = replace_na(exp_num, "not_applicable")) %>% 
  mutate(day = replace_na(day, "not_applicable")) %>% 
  # Make sure variables of interest are treated as factors
  mutate(unique_mouse_id = factor(unique_mouse_id, levels = mouse_id_levels),
         unique_cage_no = factor(unique_cage_no, levels = cage_levels),
         group = factor(group, levels = group_levels),
         exp_num = factor(exp_num, levels = exp_levels),
         day = factor(day, levels = day_levels),
         sample_type = factor(sample_type, levels = c("stool", "cecum", "proximal_colon", "distal_colon", "FMT", "water"))) #include the additional factor levels added above

#PERMANOVA (adonis) of all relevant variables
all_adonis <- adonis(all_dist~(group/(sample_type*exp_num*unique_cage_no*ext_plate*miseq_run))*day, strata = all_variables$unique_mouse_id, data = all_variables, permutations = 10)
#How to account for sample_type?
#Error: vector memory exhausted (limit reached?)

#Check subset of variables
#Test for sample_type, ext_plate, and any interactions between them.
sample_miseq_ext_plate_adonis <- adonis(all_dist~ sample_type*(miseq_run/ext_plate), data = all_variables, permutations = 1000)
#Write PERMANOVA (adonis) results to .tsv
tibble(effects = c("sample_type", "miseq_run", "miseq_run:ext_plate", "sample_type:miseq_run:ext_plate", "Residuals"),
                             r_sq = sample_miseq_ext_plate_adonis$aov.tab$R2[1:5],
                             p = sample_miseq_ext_plate_adonis$aov.tab$Pr[1:5]) %>% 
  write_tsv("exploratory/notebook/adonis_sample_miseq_plate.tsv") 

#Plot PCoA data with different variables of interest assigned to the color aesthetic----
plot_pcoa_variable <- function(color_variable){
  all_pcoa_data %>%   
    left_join(prep_variables, by = "unique_label") %>% #Add variables for ext plate number and miseq run
    ggplot(aes(x=axis1, y=axis2, color = {{ color_variable }}))+
    geom_point(size=2)+
    theme_classic()+
    theme(legend.position = "bottom")+
    xlim(-0.425, 0.65)+
    ylim(-0.525, 0.5)+
    labs(x = paste("PCoA 1 (", all_axis1, "%)", sep = ""), #Annotations for each axis from loadings file
         y = paste("PCoA 2 (", all_axis2,"%)", sep = ""))
}

#PCoA plot by sample_type
plot_pcoa_variable(sample_type)+
  ggsave("exploratory/notebook/pcoa_sample_type.pdf")
plot_pcoa_variable(ext_plate)+
  ggsave("exploratory/notebook/pcoa_ext_plate.pdf")
plot_pcoa_variable(miseq_run)+
  ggsave("exploratory/notebook/pcoa_miseq_run.pdf")
