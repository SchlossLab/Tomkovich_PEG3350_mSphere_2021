Download the [latest release](https://github.com/SchlossLab/new_project/releases/latest) to the directory and decompress


## TITLE OF YOUR PAPER GOES HERE

YOUR PAPER'S ABSTRACT GOES HERE




### Overview

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- study.Rmd    # executable Rmarkdown for this study, if applicable
	| |- study.md     # Markdown (GitHub) version of the *.Rmd file
	| |- study.tex    # TeX version of *.Rmd file
	| |- study.pdf    # PDF version of *.Rmd file
	| |- header.tex   # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- XXXX.csl     # csl file to format references for journal XXX
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
* mothur (v1.XX.0) should be located in the user's PATH
* R (v. 3.X.X) should be located in the user's PATH
* R packages:
  * `knitr`
  * `rmarkdown`
* etc


#### Running analysis
Download 16S rRNA sequencing dataset from the NCBI Sequence Read Archive (BioProject Accession no. _______).
```
git clone ________
```
Transfer 16S rRNA sequencing fastq.gz files into  ________
```
cd  ________
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
```
```
make write.paper
```
