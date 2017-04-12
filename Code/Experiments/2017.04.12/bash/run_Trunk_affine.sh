#!/bin/bash

#SBATCH
#SBATCH --job-name=run_Trunk_affine
#SBATCH --array=0-6
#SBATCH --time=2-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=120000
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

CLASSIFIERS=(rf rerf rerfr frc frcr rr_rf rr_rfr)

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
C=${CLASSIFIERS[${SLURM_ARRAY_TASK_ID}]}
matlab -nosplash -nodisplay -singleCompThread -r "run_Trunk_affine_vary_n_${C};exit" > run_Trunk_affine_vary_n_${C}.out
