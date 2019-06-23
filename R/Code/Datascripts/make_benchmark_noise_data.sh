#!/bin/bash
#SBATCH --job-name=make_benchmark_noise_data
#SBATCH --time=0:10:00
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
#---------------------------------------------------------------------
# SLURM job script to run serial R
#---------------------------------------------------------------------

ml R
ml # confirm modules used
Rscript make_benchmark_noise_data.R
