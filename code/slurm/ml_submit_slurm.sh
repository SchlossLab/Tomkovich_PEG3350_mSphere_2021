#!/bin/bash

#SBATCH --job-name=mikropml_d10_status_w_d5_genus_data_minus_WMR

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100MB
#SBATCH --time=48:00:00

#SBATCH --output=log/hpc/slurm-%j_%x.out

#SBATCH --account=pschloss1
#SBATCH --partition=standard

#SBATCH --mail-user=tomkoset@umich.edu
#SBATCH --mail-type=BEGIN,END

# To run: conda activate smk-ML before submitting job.

time snakemake --profile config/slurm --latency-wait 90 --configfile config/config.yml
