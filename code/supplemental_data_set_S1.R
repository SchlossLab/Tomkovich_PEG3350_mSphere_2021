source("code/utilities.R")

#Generate single excel workbook containing sheets for Data Set S1

#Sheet 1: 5-day PEG stool PERMANOVA results----
permanova_5_d_stools <- read_tsv("data/process/5_days_PEG_permanova_stools.tsv")


#Combine all the above tables as sheets---
sheets_combined <- list("Sheet 1" = permanova_5_d_stools, 
                    "Sheet 2" = _____,
                    "Sheet 3" = _____,
                    "Sheet 4"  = _____, 
                    "Sheet 5" = _____, 
                    "Sheet 6" = _____,
                    "Sheet 7" = _____,
                    "Sheet 8" = _____, 
                    "Sheet 9" = _____,
                    "Sheet 10" = _____,
                    "Sheet 11" = _____,
                    "Sheet 12" = _____, 
                    "Sheet 13" = _____,
                    "Sheet 14" = _____,
                    "Sheet 15" = _____)

#Write sheets to  Data Set S1 excel file----
write_xlsx(sheets_combined, "submission/Data Set S1.xlsx", format_headers = TRUE)
