library(tidyverse)

seq_count_data <- read_tsv("data/mothur/peg3350.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.pick.count.summary", col_names=c("sample", "nseqs")) %>% 
  filter(!str_detect(sample, "water")) #Removes water control samples
  
#If we rarefy to 1000 sequences per sample:
n_1000 <- seq_count_data %>% filter(nseqs < 1000) %>% select(sample) %>% nrow()
# We'll lose 0 samples

#If we rarefy to 1500 sequences per sample
n_1500 <- seq_count_data %>% filter(nseqs < 1500) %>% select(sample) %>% nrow()
# We'll lose 1 sample

#If we rarefy to 2000 sequences per sample
n_2000 <- seq_count_data %>% filter(nseqs < 2000) %>% select(sample) %>% nrow()
# We'll lose 3 samples

#Choose 2000 sequences per sample as current cutoff

#If we rarefy to 3000 sequences per sample
n_3000 <- seq_count_data %>% filter(nseqs < 3000) %>% select(sample) %>% nrow()
# We'll lose 10 samples
