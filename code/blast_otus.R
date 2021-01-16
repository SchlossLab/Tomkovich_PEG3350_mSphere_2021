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
