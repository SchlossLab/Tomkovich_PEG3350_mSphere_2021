source("code/utilities.R") #Loads libraries, reads in metadata, functions
#Identify the bacterial OTUs that correspond to C. difficile

#Read in taxonomy:
taxonomy <- read_tsv(file="data/process/peg3350.taxonomy") %>% 
  select(-Size) %>%
  mutate(Taxonomy=str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>%
  mutate(Taxonomy=str_replace_all(Taxonomy, ";$", "")) %>%
  separate(Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';')

#Sequences for Potential C. difficile OTUs:
c_diff_otus <- taxonomy %>% 
  filter(genus == "Peptostreptococcaceae_unclassified") %>% 
  pull(OTU) 
#List of 30 OTUs
#Otu0012, shows up most frequently in our dataset

#File containing the representative sequences for the OTUs identified in the dataset
#Function to read in fasta files, based on library(devtools) source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")
#I modified the function to remove the lines that export as .csv
FastaToTabular <- function (filename){
  
  #read fasta file
  
  file1 <- readLines(filename)
  
  #find the genename location by grepping >
  
  location <- which((str_sub(file1,1,1))==">")
  
  #start an empty vector to collect name and sequence 
  
  name=c()
  sequence =c()
  
  
  
  #number of genes= number of loops
  #extract name first
  for ( i in 1:length(location)){
    name_line = location[i]
    name1 = file1[name_line]
    name=c(name,name1)
    #extract sequence between the names
    #the last sequence will be missed using this strategy 
    #so, we are using if condition to extract last sequence 
    start= location[i]+1
    end = location[i+1]-1
    if ( i < length (location)){
      
      end=end
      
    } else {
      
      end=length(file1)
    }
    
    lines = start:end
    sequence1= as.character(paste(file1[lines],collapse = ""))
    sequence =c(sequence,sequence1)
  }
  
  #now create table using name and sequence vector 
  
  data <- tibble(name,sequence)
  
  
  #function ends
}

#Read in fasta file containing the representative sequences for each OTU in dataset
otu_seqs <- FastaToTabular("data/mothur/peg3350.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.rep.fasta")


#Subset otu_seqs to just the rows that have a potential c_diff otu in their name (c_diff_otus)
c_diff_seq_all <- map_df(c_diff_otus, function(c_diff_otus){
  c_diff_seq <- otu_seqs %>% 
    filter(str_detect(name, c_diff_otus)) 
})

c_diff_seq_all %>% pull(sequence)
c_diff_seq_all <- c_diff_seq_all %>% 
  mutate(ncbi_blast_result = "")

#Blast search for 1 OTU: 
#Select standard nucleotide BLAST
#For database: select rRNA/ITS databases and 16S ribosomal RNA sequences

#BLAST search for 530 Peptostreptococcaceae OTUs against C. difficile rRNA:
#Select align 2 or more sequences. Paste 30 OTU sequences in top box
#Place C. difficile rRNA gene in bottom box (used C. difficile ATCC 9689 Accession # NR_113132.1)
#Saved results in "data/process/59OTus_vs_C.diff_ATCC9689-Alignment-HitTable.csv"
blast_results <- read_csv("data/process/59OTus_vs_C.diff_ATCC9689-Alignment-HitTable.csv", 
                          col_names = c("query_acc.ver", "subject_acc.ver", "%identity", "alignment", "length", "mismatches",
                                        "gap opens", "q.start", "q.end", "subject", "evalue", "bit score")) %>%
  mutate(otu_list_no = 1:30)
#e value = Expect value parameter. Number of hits one can expect to see by chance
#bit score: sequence similarity independent of query sequence length and database size. Normalized based on the rawpairwise alignment score
#C.diff OTU list:
c_diff_otu_list <- as.data.frame(c_diff_otus) %>% 
  mutate(otu_list_no = 1:30) 

percent_identity_dist <- blast_results %>% 
  left_join(c_diff_otu_list, by = "otu_list_no") %>% 
  mutate(c_diff_otus = str_replace_all(c_diff_otus, "Otu", ""),
         c_diff_otus = str_remove(c_diff_otus, "^0+")) %>% 
  ggplot(aes(x=query_acc.ver, y = `%identity`, color = c_diff_otus, shape=c_diff_otus, show.legend = FALSE))+
  geom_text(aes(label = c_diff_otus), position = position_jitter(width = 0.5, height = 0.5))+
  scale_shape_identity()+
  labs(title="Blastn to C. difficile 16S rRNA", 
       x=NULL,
       y="% Identity")+
  theme_classic()+
  theme(plot.title = element_text(hjust =0.5),
        legend.position = "none", #Remove legend
        axis.text.x = element_blank())
save_plot("results/figures/otus_peptostreptococcaceae_blast_results.png", percent_identity_dist, base_height =5, base_width = 6)


#Important genera for FMT success as discussed in Kazemian et. al (2020) ----
#Check if the Oscillibacter genera match the species Oscillibacter sp. PEA192 described in paper's supplementary data
oscillibacter_otus <- taxonomy %>%
  filter(genus == "Oscillibacter") %>%
  pull(OTU)

oscillibacter_seq_all <- map_df(oscillibacter_otus, function(oscillibacter_otus){
  oscillibacter_seq <- otu_seqs %>% 
    filter(str_detect(name, oscillibacter_otus)) 
}) 

oscillibacter_seq_all %>% pull(sequence)
oscillibacter_seq_all <- oscillibacter_seq_all %>% 
  mutate(ncbi_blast_result = "")
#Compare the 41 OTU sequences to Oscillibacter sp. PEA192 (Accession # NZ_AP018532.1)
#Follow same parameters as described in C. diff comparison above
#Bring in BLAST results in csv format
oscillibacter_blast_results <- read_csv("data/process/41OTUs_vs_FMT_Successful_Osciliibacter.csv",
                                        col_names = c("query_acc.ver", "subject_acc.ver", "%identity", "alignment", "length", "mismatches",
                                                      "gap opens", "q.start", "q.end", "subject", "evalue", "bit score")) %>%
  mutate(otu_list_no = 1:123)

#Transform vector to dataframe for plot
oscillibacter_otu_list <- as.data.frame(oscillibacter_otus) %>%
  mutate(otu_list_no = 1:41)

#Plot percent identity of all Oscillibacter OTUs
oscillibacter_percent_identity_dist <- oscillibacter_blast_results %>%
  filter(otu_list_no == 1:41) %>%
  left_join(oscillibacter_otu_list, by = "otu_list_no") %>% 
  mutate(oscillibacter_otus = str_replace_all(oscillibacter_otus, "Otu", ""),
         oscillibacter_otus = str_remove(oscillibacter_otus, "^0+")) %>% 
  ggplot(aes(x=query_acc.ver, y = `%identity`, color = oscillibacter_otus, shape=oscillibacter_otus, show.legend = FALSE))+
  geom_text(aes(label = oscillibacter_otus), position = position_jitter(width = 0.5, height = 0.5))+
  scale_shape_identity()+
  labs(title="Blastn to Oscillibacter sp. PEA19 16S rRNA", 
       x=NULL,
       y="% Identity")+
  theme_classic()+
  theme(plot.title = element_text(hjust =0.5),
        legend.position = "none", #Remove legend
        axis.text.x = element_blank())
save_plot("results/figures/otus_oscillibacter_blast_results.png", oscillibacter_percent_identity_dist, base_height =5, base_width = 6)


#BLAST OTU 12 to check top matches----
#Pull OTU 12 from taxonomic data file
otu_12 <- taxonomy %>%
  filter(OTU == "Otu0012") %>%
  pull(OTU)
#Pull corresponding sequence of OTU 12
otu_12_seq <- map_df(otu_12, function(otu_12) {
  otu12_seq <- otu_seqs %>% 
    filter(str_detect(name, otu_12)) 
})
otu_12_seq %>% pull(sequence)
#Use single query search on otu_12 sequence
#Paste sequence in top box, selecting 16s rRNA/ITS database
#Bring in BLAST hit table results in csv format  
otu_12_blast_results <- read_csv("data/process/OTU_12_HitTable.csv")

top_otu_12_hits <- otu_12_blast_results %>% 
  filter(`Per. ident` > 96.0) %>%
  ggplot(aes(x=`Scientific Name`, y = `Per. ident` , color = `Scientific Name`, shape=`Scientific Name`, show.legend = FALSE))+
  geom_text(aes(label = `Scientific Name`), position = position_jitter(width = 0.5, height = 0.5))+
  scale_shape_identity()+
  labs(title="Top OTU 12 Blastn Results", 
       x=NULL,
       y="% Identity")+
  theme_classic()+
  theme(plot.title = element_text(hjust =0.5),
        legend.position = "none", #Remove legend
        axis.text.x = element_blank())
save_plot("results/figures/otu_12_top_blast_results.png", top_otu_12_hits, base_height =5, base_width = 6)
#BLAST all porphyromonadaceae OTUs against S24-7 bacteria (Muriibaculaceae) Accession Number: CP015402----
porphyromonadaceae_otus <- taxonomy %>% 
  filter(genus == "Porphyromonadaceae_unclassified") %>%
  pull(OTU)

#Subset otu sequences to those that have porphyromonadaceae OTUs in them
porphyromonadaceae_seq_all <- map_df(porphyromonadaceae_otus, function(porphyromonadaceae_otus) {
  porphyromonadaceae_seq <- otu_seqs %>%
    filter(str_detect(name, porphyromonadaceae_otus))
})

#Copy sequences to clipboard for pasting into BLAST
zz <- pipe('pbcopy', 'w')
porphyromonadaceae_seq_all$sequence %>%
  write.table(file = zz, sep='\t', row.names = FALSE)
close(zz)

#Select align two or mmore sequences
#Paste sequences in top box and CP015402 in bottom box 
#Bring in BLAST hit table csv
porphy_blast_results <- read_csv("data/process/porphyromonadaceae_OTUs_hitTable.csv", 
                                 col_names = c("query_acc.ver", "subject_acc.ver", "%identity", "alignment", "length", "mismatches",
                                               "gap opens", "q.start", "q.end", "subject", "evalue", "bit score")) %>%
  mutate(otu_list_no = 1:5300)
#e value = Expect value parameter. Number of hits one can expect to see by chance
#bit score: sequence similarity independent of query sequence length and database size. Normalized based on the rawpairwise alignment score
#Porphyromonadaceae OTU list
porphyr_otu_list <- as.data.frame(porphyromonadaceae_otus) %>%
  mutate(otu_list_no = 1:1329)

#Plot percent identity of all otus over 97%, plot too cluttered with all 1329 otus
porphyr_percent_identity_dist <- porphy_blast_results %>% 
  filter(`%identity` > 97.0) %>%
  left_join(porphyr_otu_list, by = "otu_list_no") %>% 
  mutate(porphyromonadaceae_otus = str_replace_all(porphyromonadaceae_otus, "Otu", ""),
         porphyromonadaceae_otus = str_remove(porphyromonadaceae_otus, "^0+")) %>% 
  ggplot(aes(x=query_acc.ver, y = `%identity`, color = porphyromonadaceae_otus, shape=porphyromonadaceae_otus, show.legend = FALSE))+
  geom_text(aes(label = porphyromonadaceae_otus), position = position_jitter(width = 0.5, height = 0.5))+
  scale_shape_identity()+
  labs(title="Blastn to M. intestinale strain YL27 16S rRNA", 
       x=NULL,
       y="% Identity")+
  theme_classic()+
  theme(plot.title = element_text(hjust =0.5),
        legend.position = "none", #Remove legend
        axis.text.x = element_blank())
save_plot("results/figures/otus_porphyromonadaceae_blast_results.png", porphyr_percent_identity_dist, base_height =5, base_width = 6)


#BLAST all bacteriodales OTUs against S24-7 bacteria (Muriibaculaceae) Accession Number: CP015402----
#Pull bacteriodales OTUs
bacteriodales_otus <- taxonomy %>%
  filter(genus == "Bacteroidales_unclassified") %>%
  pull(OTU)
  
#Subset otu sequences with matches to the OTUs in bacteriodales otu list
bacteriodales_seq_alll <- map_df(bacteriodales_otus, function(bacteriodales_otus) {
  bacteriodales_seq <- otu_seqs %>%
    filter(str_detect(name, bacteriodales_otus))
})
#Create pipe to copy all 126 sequences to clipboard for pasting into BLAST
zz <- pipe('pbcopy', 'w')
bacteriodales_seq_alll$sequence %>%
  write.table(file = zz, sep='\t', row.names = FALSE)
close(zz)

#Select align two or mmore sequences
#Paste sequences in top box and CP015402 in bottom box 
#Bring in BLAST hit table csv
bacteriodales_blast_results <- read_csv("data/process/bacteriodales_otus_HitTable.csv",
                                        col_names = c("query_acc.ver", "subject_acc.ver", "%identity", "alignment", "length", "mismatches",
                                                      "gap opens", "q.start", "q.end", "subject", "evalue", "bit score")) %>%
  mutate(otu_list_no = 1:480)

bacteriodales_otu_list <- as.data.frame(bacteriodales_otus) %>%
  mutate(otu_list_no = 1:126)

bacteriodales_percent_identity_dist <- bacteriodales_blast_results %>% 
  filter(otu_list_no == 1:126) %>%
  left_join(bacteriodales_otu_list, by = "otu_list_no") %>% 
  mutate(bacteriodales_otus = str_replace_all(bacteriodales_otus, "Otu", ""),
         bacteriodales_otus = str_remove(bacteriodales_otus, "^0+")) %>% 
  ggplot(aes(x=query_acc.ver, y = `%identity`, color = bacteriodales_otus, shape=bacteriodales_otus, show.legend = FALSE))+
  geom_text(aes(label = bacteriodales_otus), position = position_jitter(width = 0.5, height = 0.5))+
  scale_shape_identity() +
  labs(title="Blastn to M. intestinale strain YL27 16S rRNA", 
       x=NULL,
       y="% Identity")+
  theme_classic()+
  theme(plot.title = element_text(hjust =0.5),
        legend.position = "none", #Remove legend
        axis.text.x = element_blank())
save_plot("results/figures/otus_bacteriodales_blast_results.png", bacteriodales_percent_identity_dist, base_height =5, base_width = 6)

