#!/bin/bash

#SBATCH
#SBATCH --job-name=run_orthant
#SBATCH --time=7-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=120000
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

matlab -nosplash -nodisplay -singleCompThread -r "run_Orthant_rerf;exit"

echo "job complete"
