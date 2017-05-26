#!/bin/bash

#SBATCH
#SBATCH --job-name=run_benchmarks
#SBATCH --array=1-106
#SBATCH --time=2-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=120000
#SBATCH --partition=bw-parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

NAME_FILE=~/work/tyler/Benchmarks/Data/uci/Names.txt
DATASET=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $NAME_FILE)

matlab -nosplash -nodisplay -singleCompThread -r "run_${DATASET}_2017_05_27;exit"

echo "job complete"
