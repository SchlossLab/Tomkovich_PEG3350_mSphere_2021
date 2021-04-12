source("code/utilities.R") #Loads libraries, reads in metadata, functions

library(devtools)
#install_github("zdk123/SpiecEasi") how to install SpiecEasi
#Initially: got the following error:
#clang: error: linker command failed with exit code 1 (use -v to see invocation)
#make: *** [SpiecEasi.so] Error 1
#ERROR: compilation failed for package ‘SpiecEasi’
#* removing ‘/Library/Frameworks/R.framework/Versions/4.0/Resources/library/SpiecEasi’
#Error: Failed to install 'SpiecEasi' from GitHub: 
#Solution:
#Similar error encountered by others with potential solutions: https://github.com/zdk123/SpiecEasi/issues/138
#Had to install fortran: https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/
#Alternative option: install package with conda: https://anaconda.org/bioconda/r-spieceasi 
library(SpiecEasi)
library(Matrix)
library(igraph)
library(ggnet) #To install: devtools::install_github("briatte/ggnet")
library(GGally)
library(network)
library(sna)
source('code/functions/ggnet2.R') #Modified by Nick in his 2020 preprint paper
#Nick's implemnetation of Spiec-Easi: https://github.com/SchlossLab/Lesniak_Clearance_XXXX_2020/blob/master/code/build_fig6.R

set.seed(19760620) #Same seed used for mothur analysis
seed <- 19760620 

metadata <- metadata %>%
  mutate(day = as.integer(day))  #Day variable (transformed to integer to get rid of decimals on PCoA animation

#Analysis of project data
#Use subsample shared file (removes samples with < 1000 sequences)
otu_data <- read_tsv("data/process/peg3350.opti_mcc.0.03.subsample.shared") %>%
  select(-label, -numOtus) %>%
  rename(unique_label = Group) %>% 
  left_join(select(metadata, group, unique_label, day, avg_cfu), by = "unique_label") #Join to metadata

#Read in taxonomy labels
network_labels <- read_tsv(file="data/process/peg3350.taxonomy") %>%
  select(-Size) %>%
  mutate(key=OTU) %>%
  mutate(key=str_to_upper(key)) %>%
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>%
  mutate(taxa=gsub("(.*)_.*","\\1",Taxonomy)) %>%
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>%
  mutate(taxa=str_replace_all(taxa, c("Clostridium_" = "Clostridium "))) %>% 
  mutate(taxa=gsub(".*;","",taxa)) %>%
  mutate(taxa=gsub("(.*)_.*","\\1",taxa)) %>%
  mutate(taxa=gsub('[0-9]+', '', taxa)) %>%
  mutate(taxa=str_remove_all(taxa, "[(100)]")) %>%
  unite(key, taxa, key, sep=" (") %>%
  mutate(key = paste(key,")", sep="")) %>%
  select(-OTU, -Taxonomy) %>%
  rename(otu=key) %>%
  mutate(otu=paste0(gsub('TU0*', 'TU ', otu))) %>%
  separate(otu, into = c("bactname", "OTUnumber"), sep = "\\ [(]", remove = FALSE) %>% #Add columns to separate bacteria name from OTU number to utilize ggtext so that only bacteria name is italicized
  mutate(otu_name = glue("*{bactname}* ({OTUnumber}")) %>% #Markdown notation so that only bacteria name is italicized
  pull(otu_name)

#Arguments for spiec-easi
se_pargs <- list(rep.num=99, seed = seed, ncores = 4)
# function to set offset angle to arrange node labels outside circle
#  kjhealy/polar-labels.r https://gist.github.com/kjhealy/834774
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

#Test of spiec-easi
se_peg <- otu_data %>% 
  filter(group %in% c("WM", "WMC", "WMR",
                      "M1", "1RM1",
                      "CWM", "FRM", "RM") & day > 0) %>% 
  select(-unique_label, -group, -day, Cdiff = avg_cfu) %>% 
  filter(!is.na(Cdiff))
se_peg <- otu_data %>% 
  filter(group == "WM" & day > 0) %>% 
  select(-unique_label, -group, -day, Cdiff = avg_cfu) %>% 
  filter(!is.na(Cdiff))

otus_present <- se_peg %>% 
  summarise_all(function(x){
                  sum(x >= 1) >= .1 * nrow(se_wm) #select OTUS present in 10% of samples for the group of interest
          }) %>% 
  gather(OTU, present) %>% 
  filter(present == TRUE) %>% 
  pull(OTU)

#Taking two long with all OTUs. Limit to 10% of OTUs in the selected group (69 OTUs total)
se_peg <- se_peg %>% 
  select(otus_present) %>% 
  as.matrix

se_model <- spiec.easi(se_peg, method = 'mb', lambda.min.ratio = 1e-3, nlambda = 500,
                       sel.criterion = 'bstars', pulsar.select = TRUE, pulsar.params = se_pargs)
se_df <- se_peg
#create igraph object
ig.mb <- adj2igraph(getRefit(se_model))
# set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(se_df, 1))+6
se_coord <- layout.fruchterman.reingold(ig.mb)
par(mfrow=c(1,1))
wm_plot <- plot(ig.mb, layout = se_coord, vertex.size = vsize, vertex.label=NA, main = "WM")

# determine edge weights
se_beta <- symBeta(getOptBeta(se_model), mode='maxabs')
elist.mb <- summary(sebeta)
hist(elist.mb[,3], main='', xlab = 'edge weights')

#Degree statistics from the network
dd.mb <- degree.distribution(ig.mb)
plot(0:(length(dd.mb)-1), dd.mb, ylim=c(0,.35), type='b',
     ylab="Frequency", xlab="Degree", main="Degree Distributions")

# set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(se_df, 1))+6
#Determine edge weights
se_beta <- symBeta(getOptBeta(se_model), mode='maxabs')
se_edges <- Matrix::summary(se_beta)
# network degree distributions
se_network <- adj2igraph(getRefit(se_model))
se_dd <- degree.distribution(se_network)
# determine network stability (closest to 0.05 is best, increase nlambda if not close)
se_stability <- getStability(se_model)
if(se_stability < 0.045){ stop(paste0('Stability low (', se_stability, 
                                      '), increase nlambda or decrease lambda.min.ratio arguments'))}

## setup model output for graphing network
se_interaction_matrix <- getRefit(se_model)
# name matrix positions
colnames(se_interaction_matrix) <- rownames(se_interaction_matrix) <- gsub('Otu0*', '', colnames(se_df))
# subset network to only OTUs directly interacting with C. difficile
first_order_otus <- c(names(which(se_interaction_matrix[,'Cdiff'] > 0)), 'Cdiff')
#second_order_otus <- names(apply(se_data[,otus], 1 , sum) > 0)
cdiff_interactions <- se_interaction_matrix[first_order_otus, first_order_otus]
labels <- c(network_labels[as.numeric(head(colnames(cdiff_interactions), -1))], '*C. difficile*')
colnames(cdiff_interactions) <- rownames(cdiff_interactions) <- labels
#se_interaction_matrix <- se_interaction_matrix[second_order_otus, second_order_otus]
# add edge weights
# names matrix positions
colnames(se_beta) <- rownames(se_beta) <- gsub('Otu0*', '', colnames(se_df))
# subset network to OTUs interacting with C difficile
wt_cdiff_interactions <- se_beta[first_order_otus, first_order_otus]
colnames(wt_cdiff_interactions) <- rownames(wt_cdiff_interactions) <- labels
# create igraph object with edge weights
wt_first_order_network <- adj2igraph(wt_cdiff_interactions, 
                                     vertex.attr = list(name = colnames(wt_cdiff_interactions)))

# setup network attributes to create igraph network graph for output
vsize <- vsize[gsub('Otu0*', '', names(vsize)) %in% first_order_otus]
edge_wt <- abs(E(wt_first_order_network)$weight)
edge_direction <- ifelse(E(wt_first_order_network)$weight < 0, 'red', 'blue')
lab.locs <- radian.rescale(x=1:length(first_order_otus), direction=-1, start=0)
network <- adj2igraph(cdiff_interactions, 
                      vertex.attr = list(name = colnames(cdiff_interactions), size = vsize^2/10, 
                                         color = "purple", label.color='black', label.cex = 0.7, label.dist = 2, 
                                         label.degree = lab.locs),
                      edge.attr = list(width = edge_wt*10, color = edge_direction))

se_output <- list(edge_wts = se_edges, degree_dist = se_dd, stability = se_stability, 
                  all_otus = colnames(se_df), otus = first_order_otus, full_network = se_network, 
                  cdiff_network = network)
return(se_output)

#Simple plot of network
plot(network)
plot(se_output$cdiff_network)
class(network) #object is class igraph. 5 vertices and 6 edges
E(network) #Network has 6 edges
V(network) #Network has 5 vertices
network[] #Directly examine the network matrix
edge_attr(network) #Examine edge attributes
vertex_attr(network) #Examine vertex attributes
peg_network_graph <- ggnet2(se_output$cdiff_network, mode = 'kamadakawai',
                               color = 'color', label = T, size = 'size', vjust = 1.3, label.size = 3.5,
                               edge.size = 'width', edge.color = 'color', layout.exp = 0.2) +
  guides(size = FALSE) +
  ylim(-0.1, 1) + xlim(-0.1, 1.1) + 
  theme(axis.title = element_blank(), axis.text =element_blank(),
        axis.ticks = element_blank())


networks <- list(clinda_network, strep_network, cef_network, strep_colonized_network, cef_colonized_network)
# centrality
# all look fairly similar
# slightly lower amount of high degree for comminities remaining colonized
# Cef has significantly different betweenness, 
#  cleared communities have much higher betweenness centrality
#   so cleared communities have slightly more connections 
get_centrality <- function(x){
  net <- x$full_network
  tibble(antibiotic = x$antibiotic,
         clearance = x$clearance,
         #otu = x$all_otus,
         degree = igraph::degree(net, mode="in"), # number of its adjacent edges
         betweenness = igraph::betweenness(net, directed=T, weights=NA)) %>% # the number of shortest paths going through node
    gather(metric, value, -antibiotic, -clearance)
}

centrality_df <- map_dfr(networks, get_centrality) 
# test diffs
pvalue_df <- centrality_df %>% 
  mutate(subset = paste(antibiotic, clearance, sep = '_')) %>% 
  select(metric, subset, value) %>% 
  group_by(metric) %>% 
  nest()  %>% 
  mutate(data = map(data, function(data) nest(group_by(data, subset)))) %>% 	
  mutate(data = map(data, function(nested_df){
    test_df <- sapply(nested_df$data, function(x) sapply(nested_df$data, function(y) wilcox.test(x$value,y$value)$p.value)) %>% 
      data.frame
    test_df[upper.tri(test_df)] <- NA # set all of upper triangle to 1 to eliminate duplicate comparisons
    colnames(test_df) <- nested_df$subset
    test_df <- test_df %>% 
      mutate(row_names = nested_df$subset) %>% 
      pivot_longer(col = -row_names, names_to = 'col_names', values_to = 'pvalue') %>% 
      filter(pvalue != 1, !is.na(pvalue)) %>% # eliminate all self comparisons and upper triangle
      mutate(pvalue = p.adjust(pvalue, method = 'BH')) # correct p values
    return(test_df)
  })) %>% 
  unnest(data)

annotation_df <- data.frame(metric = rep(c('Degree Centrality', 'Betweenness Centrality'), each = 6),
                            x1=c(1, 1, 2.85, 2.85, 3.15, 3.15, 1, 1, 2.85, 2.85, 3.15, 3.15), 
                            x2=c(1.85, 2.15, 1.85, 2.15, 1.85, 2.15, 1.85, 2.15, 1.85, 2.15, 1.85, 2.15), 
                            xnote = c(1.5, 1.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 2.5, 2.5, 2.5, 2.5),
                            y1 = c(10^0.843, 10^0.878, 10^0.913, 10^0.948, 10^0.983, 10^1.018, 
                                   10^2.41, 10^2.51, 10^2.61, 10^2.71, 10^2.81, 10^2.91),
                            ynote = c(10^0.8445, 10^0.8795, 10^0.9145, 10^0.9495, 10^0.9845, 10^1.0195, 
                                      10^2.42, 10^2.52, 10^2.62, 10^2.72, 10^2.82, 10^2.92),
                            annotations='*', clearance = 'NA', antibiotic = 'NA')

centrality_plot <- centrality_df %>% 
  mutate(antibiotic = factor(antibiotic, levels = c('Clindamycin', 'Cefoperazone', 'Streptomycin')),
         metric = case_when(metric == 'betweenness' ~ 'Betweenness Centrality',
                            metric == 'degree' ~ 'Degree Centrality',
                            T ~ metric)) %>% 
  ggplot(aes(x = antibiotic, y = value + .1, fill = clearance, color = antibiotic)) + 
  geom_boxplot(width = 0.6,  position = position_dodge2(preserve = "single"),
               show.legend = F) + 
  geom_point(size = NA, shape = 22) + 
  scale_color_manual(values = abx_color$color, limits = abx_color$abx) + 
  scale_fill_manual(values = c(NA, 'gray'), limits = c('Cleared', 'Colonized')) + 
  facet_wrap(.~metric, scales = 'free') + 
  scale_y_log10(breaks = c(0.1, 1, 10, 100),
                labels = c('0', '1', '10', '100')) + 
  guides(color = 'none') + theme_bw() + 
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(legend.position = c(0.5, 0.565),
        panel.spacing = unit(5,'lines'),
        strip.background = element_blank(),
        strip.text = element_text(size = 14)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) + 
  geom_text(data = annotation_df, aes(x = xnote, y = ynote, label = annotations), 
            color = 'black') +
  geom_segment(data = annotation_df, aes(x = x1, xend = x2, y = y1, yend = y1), 
               color = 'black', size = 0.25)


#igraph tutorial---- https://kateto.net/netscix2016.html 
gl <- graph(edges=c(1,2, 2,3, 3,1), n =3, directed =F)
plot(gl)  

