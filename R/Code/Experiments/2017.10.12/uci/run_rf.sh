#!/bin/bash

#SBATCH
#SBATCH --job-name=run_uci
#SBATCH --time=4-04:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=120G
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

Rscript run_Benchmarks_rf_2017_10_12.R

echo "job complete"
