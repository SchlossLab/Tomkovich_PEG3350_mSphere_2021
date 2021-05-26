source("code/utilities.R")

#Generate single excel workbook containing sheets for Data Set S1

#Sheet 1: 5-day PEG subset stool PERMANOVA results
permanova_5_d_stools <- read_tsv("data/process/5_days_PEG_permanova_stools.tsv")

#Sheet 2: 5-day PEG subset stool Shannon diversity analysis
f_d_group_shannon <- read_tsv("data/process/5_days_PEG_shannon_stats_stools_subset.tsv")

#Sheet 3: 5-day PEG impacted genera in stool communities
PEG_genera <- read_tsv("data/process/5_days_PEG_genus_PEG_paired.tsv")

#Sheet 4: Genera with different relative abundances between 5-day PEG subset treatment groups in stool communities
f_d_group_genera <- read_tsv("data/process/5_days_PEG_genus_group_stools.tsv")

#Sheet 5: Statistical test of genera in WMR mice between day 1 & day 8 time points
WMR_genera <- read_tsv("data/process/5_days_PEG_genus_WMR_paired.tsv")

#Sheet 6: 5-day PEG subset cecum Shannon diversity analysis
f_d_cecum_group_shannon <- read_tsv("data/process/5_days_PEG_shannon_stats_cecum_subset.tsv")

#Sheet 7: 5-day PEG subset tissue PERMANOVA results
permanova_5_d_tissues <- read_tsv("data/process/5_days_PEG_permanova_tissues.tsv")

#Sheets 8-10: Genera with different relative abundances between 5-day PEG subset treatment groups in tissue communities
#Sheet 8: Cecum
PEG_genera_cecum <- read_tsv("data/process/5_days_PEG_genus_group_cecum.tsv")
#Unclassified was also different besides the 4 genera mentioned in the text
#Sheet 9: Proximal colon
PEG_genera_pc <- read_tsv("data/process/5_days_PEG_genus_group_proximal_colon.tsv")
#Sheet 10: Distal colon
PEG_genera_dc <- read_tsv("data/process/5_days_PEG_genus_group_distal_colon.tsv")

#Sheet 11: 1-day PEG subset PERMANOVA results
permanova_1_d <- read_tsv("data/process/1_day_PEG_permanova.tsv")

#Sheet 12: 1-day subset PEG Shannon diversity statistics
shannon_1_d <- read_tsv("data/process/1_days_PEG_shannon_stats_subset.tsv")

#Sheet 13: Genera with different relative abundances between baseline and day 1 in 1-day PEG subset
b_d1_1_day_PEG_genera <- read_tsv("data/process/1_Day_PEG_genus_stats_BtoD1.tsv")

#Sheet 14: Genera with different relative abundances between baseline and day 7 in 1-day PEG subset
b_d7_1_day_PEG_genera <- read_tsv("data/process/1_Day_PEG_genus_stats_BtoD7.tsv")

#Sheet 15: post-challenge PEG subset Shannon diversity statistics
shannon_post_cdi <- read_tsv("data/process/post_CDI_PEG_shannon_stats_stools_subset.tsv")

#Sheet 16: post-challenge PEG subset richness statistics
richness_post_cdi <- read_tsv("data/process/post_CDI_PEG_richness_stats_stools_subset.tsv")

#Sheet 17: post-challenge PEG subset PERMANOVA results
permanaova_post_cdi <- read_tsv("data/process/post_CDI_PEG_permanova_stools.tsv")

#Sheet 18: Genera with different relative abundances between post-challenge subset treatment groups
post_cdi_group_genera <- read_tsv("data/process/post_CDI_PEG_genus_group_stools.tsv") 

#Sheet 19: AUROC results for 3 models tested
auc_results <-  read_csv("results/performance_results.csv")

#Combine all the above tables as sheets---
sheets_combined <- list("Sheet 1" = permanova_5_d_stools, 
                    "Sheet 2" = f_d_group_shannon,
                    "Sheet 3" = PEG_genera,
                    "Sheet 4"  = f_d_group_genera, 
                    "Sheet 5" = WMR_genera, 
                    "Sheet 6" = f_d_cecum_group_shannon,
                    "Sheet 7" = permanova_5_d_tissues,
                    "Sheet 8" = PEG_genera_cecum, 
                    "Sheet 9" = PEG_genera_pc,
                    "Sheet 10" = PEG_genera_dc,
                    "Sheet 11" = permanova_1_d,
                    "Sheet 12" = shannon_1_d, 
                    "Sheet 13" = b_d1_1_day_PEG_genera,
                    "Sheet 14" = b_d7_1_day_PEG_genera,
                    "Sheet 15" = shannon_post_cdi,
                    "Sheet 16" = richness_post_cdi,
                    "Sheet 17" = permanaova_post_cdi,
                    "Sheet 18" = post_cdi_group_genera,
                    "Sheet 19" = auc_results
                    )

#Write sheets to  Data Set S1 excel file----
write_xlsx(sheets_combined, "submission/Data Set S1.xlsx", format_headers = TRUE)
