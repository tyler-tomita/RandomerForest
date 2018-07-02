#!/bin/bash

#SBATCH
#SBATCH --job-name=ccf_sparse_parity
#SBATCH --time=2-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=120G
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com


/cm/shared/apps/MATLAB/R2018a-1/bin/matlab -nosplash -nodisplay -singleCompThread -r "run_sparse_parity_ccf_2018_07_02"

echo "job complete"
