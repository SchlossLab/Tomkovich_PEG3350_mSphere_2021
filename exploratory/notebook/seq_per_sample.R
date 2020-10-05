source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Examine the number of sequences per sample to determine subsampling parameter

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
#Lose 125 samples

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
