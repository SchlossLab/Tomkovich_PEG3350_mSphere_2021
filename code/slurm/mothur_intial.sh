#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=mothur_motilityp1-2

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=10:00:00

# Account
#SBATCH --account=pschloss99
#SBATCH --partition=standard

# Logs
#SBATCH --mail-user=tomkoset@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x-%j.out

# Environment
#SBATCH --export=ALL

# To run mothur: conda activate mothur_v1.43 before submitting job.

#####################
#                   #
#  2) Job Commands  #
#                   #
#####################


mothur code/get_error.batch
mothur code/get_good_seqs.batch
mothur code/get_shared_otus.batch
mothur code/alpha_beta.batch
