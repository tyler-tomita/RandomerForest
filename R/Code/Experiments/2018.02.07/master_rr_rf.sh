#!/bin/bash

#SBATCH
#SBATCH --job-name=rr_rf_uci
#SBATCH --array=5,6,7,34,35,36,42,43,58,62,70,85,95
#SBATCH --time=1-0:0:0
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

sed "s/abalone/${DATASET}/g" run_Benchmarks_rr_rf_2018_02_07.R > task${SLURM_ARRAY_TASK_ID}_rr_rf.R

Rscript task${SLURM_ARRAY_TASK_ID}_rr_rf.R

rm task${SLURM_ARRAY_TASK_ID}_rr_rf.R

echo "job complete"
