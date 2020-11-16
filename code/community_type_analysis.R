source("code/utilities.R") #Loads libraries, reads in metadata, functions

#Create .shared and .taxonomy files at the genus level to use in Dirichlet Multinomial Mixture analysis

#Shared file:
shared <- read.delim('data/process/peg3350.opti_mcc.0.03.subsample.shared', header=T, sep='\t') %>%
  select(-label, -numOtus) %>% 
  gather(-Group, key=OTU, value=count)

#Read in taxonomy and select genus level:
taxonomy <- read_tsv(file="data/process/peg3350.taxonomy") %>% 
  select(-Size) %>%
  mutate(Taxonomy=str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>%
  mutate(Taxonomy=str_replace_all(Taxonomy, ";$", "")) %>%
  separate(Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';') %>%
  select(OTU, "genus") %>%
  rename(taxon = genus)

unique_taxonomy <- taxonomy %>%
  select(taxon) %>%
  unique() %>%
  mutate(otu = paste0("Otu", str_pad(1:nrow(.), width=nchar(nrow(.)), pad="0")))

#Join genus level taxonomy to shared to create shared file at the genus level:
genus_shared <- inner_join(shared, taxonomy, by="OTU") %>%
  group_by(taxon, Group) %>%
  summarize(count = sum(count)) %>%
  ungroup() %>%
  inner_join(., unique_taxonomy) %>%
  select(-taxon) %>%
  spread(otu, count) %>%
  mutate(label="genus", numOtus=ncol(.)-1) %>%
  select(label, Group, numOtus, everything())
write_tsv(genus_shared, path = "data/process/peg3350.subsample.genus.shared")

select(genus_shared, -label, -numOtus) %>%
  gather(otu, count, -Group) %>%
  group_by(otu) %>%
  summarize(count=sum(count)) %>%
  inner_join(., unique_taxonomy) %>%
  rename("OTU"="otu", "Size"="count", "Taxonomy"="taxon") %>%
  write_tsv(path ="data/process/peg3350.genus.taxonomy")

#Visualize get.communitytype analysis results----
#Read in data to evaluate community type fit depending on the number of community types
dmm_fit <- read_tsv("data/process/peg3350.subsample.genus.genus.dmm.mix.fit")

laplace_plot <- dmm_fit %>% 
  ggplot()+
  geom_line(aes(x = K, y = Laplace))+
  theme_classic()
#Save results
save_plot(filename = "exploratory/notebook/motility_community_type_laplace.png", laplace_plot)

