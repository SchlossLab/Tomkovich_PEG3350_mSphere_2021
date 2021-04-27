#Load packages
library(tidyverse)
library(lubridate)
source("code/utilities.R") #Loads libraries, reads in metadata, functions

### Join sequence samples (samples we have raw sequencing data for) to metadata
### Modify metadata to create collection date column, select columns to include as NCBI SRA metadata
metadata <- peg3350.files_unique_label %>% 
  left_join(metadata) %>% 
  mutate(collection_date = str_extract(tube_label, "[0-9/]{5,}")) %>% #Extract collection date from tube label name
  mutate(collection_date = mdy(collection_date)) %>%  #Convert character to date format required by Mimarks
  mutate(collection_date = as.character(collection_date)) %>% #Turn dates back into characters so we can add missing entries to mock and water control samples
  #select columns to be included in NCBI metadata (MIMARKS peg3350.tsv file)
  select(unique_label, unique_mouse_id, group, unique_cage_no, day, exp_num, 
         sample_type, collection_date, ext_plate, miseq_run) %>% 
  mutate(sample_type = as.character(sample_type)) %>% #turn back into character vector to use with case_when
  mutate(sample_type = case_when(str_detect(unique_label, "water") ~ "water",
                                 str_detect(unique_label, "mock") ~ "mock",
                                 str_detect(unique_label, "PBS") ~ "PBS",
                                 TRUE ~ sample_type)) %>% 
  #Create additional columns required by mimarks
  mutate(description = sample_type,
         sample_title = unique_label,
         seq_method = "http://www.mothur.org/w/images/0/0c/Wet-lab_MiSeq_SOP.pdf",
         organism = case_when(sample_type %in% c("stool", "cecum", "proximal_colon", "distal_colon") ~ "mouse gut metagenome",
                              TRUE ~ "metagenome"),
         collection_date = case_when(is.na(collection_date) ~ "missing", #Mark collection_date as missing for mock & water controls
                                     TRUE ~ collection_date),
         env_biome = case_when(sample_type %in% c("stool", "cecum", "proximal_colon", "distal_colon") ~ "mouse gut",
                               TRUE ~ "not_applicable"),
         env_feature = sample_type,
         env_material = sample_type,
         geo_loc_name = "USA: Michigan: Ann Arbor",
         host = case_when(sample_type %in% c("stool", "cecum", "proximal_colon", "distal_colon") ~ "Mus musculus",
                          TRUE ~ "not_applicable"),
         lat_lon = "42.2820973 N 83.7338904 W") %>% 
  #Reorder columns to match mimarks format
  select(unique_label, description, sample_title, seq_method, organism, collection_date, env_biome, env_feature,
         env_material, geo_loc_name, host, lat_lon, unique_mouse_id, group, unique_cage_no, day, exp_num, 
         sample_type, ext_plate, miseq_run)

#Mimarks peg.tsv generated using mothur with the following command in the terminal after initiating mothur:
#get.mimarkspackage(inputdir=data/raw, file=peg3350.files, package = host_associated, requiredonly=t)

#Import MIMARKS peg3350.tsv file into data frame
mimarks <- read_tsv("data/raw/peg3350.tsv", skip = 7) 
#Add metadata columns to mimarks file to finish filling out the information for all samples
mimarks <- left_join(mimarks, metadata, by = c("#{sample name}" = "unique_label"))
write_tsv(mimarks, 'data/raw/peg3350.tsv', na = " ")

# To prepare for NCBI submission. 
#1. Open up file in excel and move additional column names to their appropriate location 
#2. Reformat date column to YYYY-mm-dd using format cells
#3. Insert seven empty rows at the top of the file
#4. Copy the following information into rows 1-7 of file:
#This is a tab-delimited file. Additional Documentation can be found at http://www.mothur.org/wiki/MIMarks_Data_Packages.
#Please fill all the required fields indicated with '*'
#Unknown or inapplicable fields can be assigned 'missing' value.
#You may add extra custom fields to this template. Make sure all the fields are separated by tabs.
#You may remove any fields not required (marked with '*'). Make sure all the fields are separated by tabs.
#You can edit this template using Microsoft Excel or any other editor. But while saving the file please make sure to save them as 'TAB-DELIMITED' TEXT FILE.
#MIMARKS.survey.host-associated.4.0

#Fill out .project file accordingly for the sequencing project.

#Once the MIMARKS and .project file have been created you are ready to run the make.sra command in terminal after navigating to project directory and intiating mothur:
# make.sra(inputdir=data/raw, file=peg3350.files, project=peg3350.project, mimark=peg3350.tsv)

#After creating the .project file, the mimarks file, and running make.sra package, you are ready to submit.
#Instructions: https://mothur.org/wiki/creating_a_new_submission/
#I prefer using the ftp instructions and transferring the files using FileZilla