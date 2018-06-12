#!/bin/bash

#SBATCH
#SBATCH --job-name=ccf_naive_uci
#SBATCH --array=1-23,25-106
#SBATCH --time=4-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=120G
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

NAME_FILE=~/work/tyler/Data/uci/processed/names.txt
DATASET=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $NAME_FILE)

sed "s/abalone/${DATASET}/g" run_abalone_ccf_naive_cat_2018_06_12.m > task_naive_${SLURM_ARRAY_TASK_ID}.m

matlab -nosplash -nodisplay -singleCompThread -r "task${SLURM_ARRAY_TASK_ID}.m"

rm task${SLURM_ARRAY_TASK_ID}.R

echo "job complete"
