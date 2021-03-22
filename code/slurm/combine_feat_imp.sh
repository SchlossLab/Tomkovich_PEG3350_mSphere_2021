#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=combine_feat_imp

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=10-00:00

# Account
#SBATCH --account=pschloss1
#SBATCH --partition=standard

# Logs
#SBATCH --mail-user=tomkoset@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x-%j.out

# Environment
#SBATCH --export=ALL


#####################
#                   #
#  2) Job Commands  #
#                   #
#####################


bash code/cat_feat_imp_csv_files.sh
