#!/bin/bash

#SBATCH
#SBATCH --job-name=ccf_naive_uci
#SBATCH --array=44
#SBATCH --time=0-2:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=120G
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

NAME_FILE=~/work/tyler/Data/uci/processed/names.txt
DATASET=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $NAME_FILE)

sed "s/abalone/${DATASET}/g" run_abalone_noise_D10_ccf_naive_cat_2019_06_03.m > task_D10_${SLURM_ARRAY_TASK_ID}.m

matlab -nosplash -nodisplay -singleCompThread -r "task_D10_${SLURM_ARRAY_TASK_ID}"

rm task_D10_${SLURM_ARRAY_TASK_ID}.m

echo "job complete"
