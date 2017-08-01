#!/bin/bash

#SBATCH
#SBATCH --job-name=run_benchmarks
#SBATCH --time=2-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=120G
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

matlab -nosplash -nodisplay -singleCompThread -r "RerF_test;exit"

echo "job complete"
