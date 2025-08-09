#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=150G
#SBATCH --job-name=genie3_job
#SBATCH -o vaishnavi/thesis/vaishhnavi_master-thesis/slurm_updated_genie/genie3_job.%j.out
#SBATCH -e vaishnavi/thesis/vaishhnavi_master-thesis/slurm_updated_genie/genie3_job.%j.err
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny



# run the R script
Rscript vaishnavi/thesis/vaishhnavi_master-thesis/genie3.R
