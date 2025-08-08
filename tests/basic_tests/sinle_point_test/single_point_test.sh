#!/bin/bash

# Allow for running conda environment (without the sourcing the conda environment does not work)
base=${CONDA_EXE%/bin/conda}
source ${base}/etc/profile.d/conda.sh
conda activate run_psi4

# Actual script
my_run_psi4.py single_point --xyz ../../example_files/h2.xyz --exc --total_charge 0 --func pbe0
