#!/bin/bash

Rscript ~/RandomerForest/R/Code/Experiments/2017.10.01/Orthant/run_Orthant_rerf_2017_10_01.R

echo "Orthant complete"

Rscript run_Trunk_rerf_2017_10_07.R

echo "Trunk complete"
