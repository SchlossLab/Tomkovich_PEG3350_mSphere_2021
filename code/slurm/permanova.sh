#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=permanova_post

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=180GB
#SBATCH --time=7-00:00:00

# Account
#SBATCH --account=pschloss1
#SBATCH --partition=standard

# Logs
#SBATCH --mail-user=tomkoset@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x-%j.out

# Environment
#SBATCH --export=ALL

# To run r & vegan adonis(): conda activate rstats before submitting job.
# --max-ppsize 500000 suggested by Kelly to work around memory issues

#####################
#                   #
#  2) Job Commands  #
#                   #
#####################

Rscript code/post_CDI_PEG_permanova.R
