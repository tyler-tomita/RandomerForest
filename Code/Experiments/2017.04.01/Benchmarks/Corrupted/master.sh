#!/bin/bash

#SBATCH
#SBATCH --job-name=run_benchmarks
#SBATCH --array=2,3,8,10,19,26-29,34,39,47-50,53,55,56,58,85,86,87
#SBATCH --time=2-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=120000
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

NAME_FILE=~/work/tyler/RandomerForest/Data/Benchmarks/Benchmark_names.txt
DATASET=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $NAME_FILE)

matlab -nosplash -nodisplay -singleCompThread -r "run_${DATASET}_corrupted;exit"

echo "job complete"
