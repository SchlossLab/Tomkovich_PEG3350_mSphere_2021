## An osmotic laxative renders mice susceptible to prolonged *Clostridioides difficile* colonization and hinders clearance

Antibiotics are a major risk factor for *Clostridioides difficile* infections (CDIs) because of their impact on the microbiota. However, non-antibiotic medications such as the ubiquitous osmotic laxative polyethylene glycol (PEG) 3350, also alter the microbiota, but whether PEG impacts CDI susceptibility and clearance is unclear. To examine how PEG impacts susceptibility, we treated C57Bl/6 mice with 5-day and 1-day doses of 15% PEG in the drinking water and then challenged the mice with *C. difficile* 630. We used clindamycin-treated mice as a control because they consistently clear *C. difficile* within 10 days post-challenge (dpc). PEG treatment alone was sufficient to render mice susceptible and 5-day PEG-treated mice remain colonized for up to 30 dpc. In contrast, 1-day PEG treated mice were transiently colonized, clearing *C. difficile* within 7 dpc. To examine how PEG treatment impacts clearance, we administered a 1-day PEG treatment to clindamycin-treated, *C. difficile*-challenged mice. Administering PEG to mice after *C. difficile* challenge prolonged colonization up to 30 dpc. When we trained a random forest model with community data from 5 dpc, we were able to predict which mice would exhibit prolonged colonization (AUROC = 0.90). Five of the top ten bacterial features important for predicting prolonged colonization had high relative abundances in the community. Examining the dynamics of these bacterial during the post-challenge period revealed patterns in the relative abundances of *Bacteroides*, *Enterobacteriaceae*, *Porphyromonadaceae*, *Lachnospiraceae*, and *Akkermansia* that were associated with prolonged *C. difficile* colonization in PEG-treated mice.

### Overview

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- manuscript.Rmd    # executable Rmarkdown for this study, if applicable
	| |- manuscript.md     # Markdown (GitHub) version of the *.Rmd file
	| |- manuscript.tex    # TeX version of *.Rmd file
	| |- manuscript.pdf    # PDF version of *.Rmd file
	| |- header.tex   # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- mSphere.csl     # csl file to format references for journal 
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| |- mothur/      # mothur processed data
	| +- process/     # cleaned data, will not be altered once created;
	|                 # will be committed to repo
	|
	|- code/          # any programmatic code
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	| +- pictures/    # diagrams, images, and other non-graph graphics
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|
	+- Makefile       # executable Makefile for this study, if applicable


### How to regenerate this repository

#### Dependencies and locations
* Gnu Make should be located in the user's PATH
* mothur (v1.43.0) should be located in the user's PATH
* FFmpeg should be located in the user's PATH
* R (v. 4.0.2) should be located in the user's PATH
* R packages:
    * broom v0.7.0
    * tidyverse_1.3.0
    * cowplot v1.0.0
    * vegan v2.5-6
    * knitr v1.29
    * rmarkdown v2.3
    * ggpubr v.0.4.0
    * gganimate v1.0.6
    * readxl v1.3.1
    * writexl v1.3
    * glue v1.4.1
    * ggtext v0.1.0
	  * magick v2.6.0
  	* here 1.0.1
* Analysis assumes the use of 10 processors

#### Running analysis
Download 16S rRNA sequencing dataset from the NCBI Sequence Read Archive (BioProject Accession no. PRJNA727293).
```
git clone https://github.com/tomkoset/Tomkovich_PEG3350_mSphere_2020
```
Transfer 16S rRNA sequencing fastq.gz files into  Tomkovich_PEG3350_mSphere_2020/data/raw
```
cd  Tomkovich_PEG3350_mSphere_2020
```
Obtain the SILVA reference alignment from version 132 described at https://mothur.org/blog/2018/SILVA-v132-reference-files/. We will use the SEED v. 132, which contain 12,083 bacterial sequences. This also contains the reference taxonomy. We will limit the databases to only include bacterial sequences.
```
wget -N https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v132.tgz
tar xvzf Silva.seed_v132.tgz silva.seed_v132.align silva.seed_v132.tax
mothur "#get.lineage(fasta=silva.seed_v132.align, taxonomy=silva.seed_v132.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v132.pick.align, processors=8)"
mv silva.seed_v132.pick.align data/references/silva.seed.align
rm Silva.seed_v132.tgz silva.seed_v132.*
#Narrow to v4 region
mothur "#pcr.seqs(fasta=data/references/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
mv data/references/silva.seed.pcr.align data/references/silva.v4.align
```
Obtain the RDP reference taxonomy. The current version is v11.5 and we use a "special" pds version of the database files, which are described at https://mothur.org/blog/2017/RDP-v16-reference_files/.
```
wget -N https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset16_022016.pds.tgz
tar xvzf Trainset16_022016.pds.tgz trainset16_022016.pds
mv trainset16_022016.pds/* data/references/
rm -rf trainset16_022016.pds
rm Trainset16_022016.pds.tgz
```
Obtain the Zymo mock community data; note that Zymo named the 5 operon of Salmonella twice instead of the 7 operon.
```
wget -N https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip
unzip ZymoBIOMICS.STD.refseq.v2.zip
rm ZymoBIOMICS.STD.refseq.v2/ssrRNAs/*itochondria_ssrRNA.fasta #V4 primers don't come close to annealing to these
cat ZymoBIOMICS.STD.refseq.v2/ssrRNAs/*fasta > zymo_temp.fasta
sed '0,/Salmonella_enterica_16S_5/{s/Salmonella_enterica_16S_5/Salmonella_enterica_16S_7/}' zymo_temp.fasta > zymo.fasta
mothur "#align.seqs(fasta=zymo.fasta, reference=data/references/silva.v4.align, processors=12)"
mv zymo.align data/references/zymo_mock.align
rm -rf zymo* ZymoBIOMICS.STD.refseq.v2* zymo_temp.fasta
```

Run the 16S rRNA sequencing data through mothur:
```
mothur code/get_error.batch
mothur code/get_good_seqs.batch
mothur code/get_shared_otus.batch
mothur code/alpha_beta.batch

#Or if using HPC:
sbatch code/slurm/mothur_process_sequences.sh
```
Copy files generated by processing sequences through mothur and place into data/process.
```
bash transfer_mothur_outputs
```
Create ordinations for the 3 sets of experiments, run PERMANOVA analysis for each set. Note: recommend using HPC to run PERMANOVA analysis.
```
Rscript code/subset_analysis.R
mothur code/subset_dist_pcoa.batch
Rscript code/5_days_PEG_16S.R
Rscript code/1_day_PEG_permanova.R
Rscript code/post_CDI_PEG_permanova.R

#Recommend using HPC to create distance matrices, ordinations, and running 5-day PEG/post-CDI PEG PERMANOVA tests. Activate conda environment specified in slurm scripts before submitting.
sbatch code/slurm/mothur_ordinations.sh
sbatch code/slurm/permanova_5d.sh
sbatch code/slurm/post_cdi.sh
```
Create genus level .shared and .taxonomy files to use for 16S sequencing analysis.
```
Rscript code/genus_level_shared_taxonomy.R
```

Examine potential *C. difficile* sequences in the dataset.
```
bash code/get_oturep.batch
Rscript code/blast_otus.R
```
Examine *C. difficile* CFU, mouse weight, histology, and 16S rRNA sequencing data for the 5-day PEG treatment set of experiments.
```
Rscript code/5_days_PEG_cfu_weight.R
Rscript code/5_days_PEG_16S.R
Rscript code/5_days_PEG_histology_scores.R
Rscript code/relative_dist_to_baseline.R
```
Examine *C. difficile* CFU, mouse weight, and 16S rRNA sequencing data for the 1-day PEG treatment set of experiments.
```
Rscript code/1_day_PEG_cfu_weight.R
Rscript code/1_day_PEG_16S.R
```
Examine *C. difficile* CFU, mouse weight, and 16S rRNA sequencing data for the post-challenge 1-day PEG treatment set of experiments.
```
Rscript code/post_CDI_PEG_cfu_weight.R
Rscript code/post_CDI_PEG_16S.R
```
Prepare genus-level input data for machine learning analysis. Remove the genus containing *C. difficile* sequences from the input data.
```
Rscript code/ml_input_data.R
```
Run mikropml pipeline with the genus-level input data using snakemake on the HPC. 
Tip: snakemake -n (Dry run). Snakemake --unlock (If you get an error that the directory is locked). The conda environment for running mikropml is in config/mikropml_snakemake.yml. Activate enviroment before submitting job.
```
#Specify which models, k_fold, and training fraction in the config/config.yml file

#To run on HPC:
sbatch code/slurm/ml_submit_slurm.sh
sbatch code/slurm/combine_feat_imp.sh
```
Visualize machine learning AUROC and feature importance results.
```
Rscript code/ml_results.R
```
Create figures for the paper.
```
Rscript code/figure_1.R
Rscript code/figure_2.R
Rscript code/figure_3.R
Rscript code/figure_4.R
Rscript code/figure_5.R
Rscript code/figure_6.R
Rscript code/figure_7.R
Rscript code/figure_8.R
Rscript code/figure_S1.R
Rscript code/figure_S2.R
Rscript code/figure_S3.R
```
Compress figures for submission with lzw compression.
```
ls submission/*.tiff | xargs sips -s format tiff -s formatOptions lzw
```
Generate Data Set S1 as an excel workbook.
```
Rscript code/supplemental_data_set_S1.R
```
Generate the paper.
```
open submission/manuscript.Rmd and knit to Word or PDF document.
```

