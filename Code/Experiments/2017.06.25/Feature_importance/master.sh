#!/bin/bash

#SBATCH
#SBATCH --job-name=feature_importance
#SBATCH --time=4-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=120000
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

matlab -nosplash -nodisplay -singleCompThread -r "feature_importance_2017_06_25;exit"

echo "job complete"
