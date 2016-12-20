#!/bin/bash -l

#SBATCH
#SBATCH --job-name=test
#SBATCH --time=1-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=parallel
#SBATCH --exclusive
#SBATCH --mail-type=end
#SBATCH --mail-user=tmtomita87@gmail.com

matlab -nosplash -nodisplay < print_numcores.m > print_numcores.out