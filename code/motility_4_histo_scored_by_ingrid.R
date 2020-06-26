library(tidyverse)
library(readxl)
library(broom)

#For Ana to get in the right directory
#setwd("~/../Box/Tomkovich_intestinal_motility_XXXX_2019/")

motility4histology <- readxl::read_excel("data/motility_4/7.17.19_histo_scored_by_ingrid.xlsx", sheet = "Scored") %>% 
  mutate(Tissue=factor(Tissue, levels=c("cecum", "colon"))) #make sure Tissue is treated as a factor

#Define color scheme----
color_scheme <- c("#238b45", "#88419d", "#f768a1", "#8c96c6") #Adapted from http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=4
color_groups <- c("C", "WM", "WMC", "WMN")
color_labels <- c( "Clind.", "5-day PEG 3350", "5-day PEG 3350 + Clind.", "5-day PEG 3350 without infection")

motility4histology %>% filter(!is.na(Tissue)) %>% 
  mutate("Mean Summary Score" = X__1,
         Tissue = ifelse(Tissue == "cecum", "Cecum", "Colon")) %>%
  ggplot(aes(Tissue, `summary score2`, color = Group))+
  geom_boxplot(lwd = 1)+
  geom_jitter(shape=19, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1), show.legend = FALSE)+
  labs(y = "Summary Score")+
  scale_color_manual(name=NULL, 
                     breaks=color_groups, 
                     labels=color_labels, 
                     values=color_scheme)+ 
  scale_y_continuous(limits = c(0, 12), breaks = c(0, 4, 8, 12))+
  theme_classic()+
  ggsave("results/figures/motility_4_histo_by_ingrid.png")

#Kruskal_wallis test for differences arcoss groups in different tissue compartments
kruskal_wallis_histo <- motility4histology %>% 
  filter(Tissue %in% c("cecum", "colon")) %>% 
  group_by(Tissue) %>% 
  do(tidy(kruskal.test(`summary score2`~factor(Group), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)

#Significant for colon (p = .0229), but not cecum (p = 0.105)
#For colon, do pairwise.wilcox.test to determine which groups of mice are significantly different from each other
m_4_colon <- motility4histology %>% 
  filter(Tissue == "colon")
pairwise_wilcox_histo <- m_4_colon %>% 
  tidy(pairwise.wilcox.test(g = factor(m_4_colon$Group), x = m_4_colon$`summary score2`, p.adjust.method = "BH"))
#Error: C stack usage  7969280 is too close to the limit