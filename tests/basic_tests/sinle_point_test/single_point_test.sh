#!/bin/bash

activate_conda_env() {
        conda_env=${1}
        if [ -z ${conda_env} ]; then
                echo called activate_conda_env without providing environment name
        fi
        base=${CONDA_EXE%/bin/conda}
        source ${base}/etc/profile.d/conda.sh
        conda activate ${conda_env}
}

conda_env=run_psi4
activate_conda_env ${conda_env}

my_run_psi4.py single_point --xyz ../../example_files/h2.xyz --exc --total_charge 0 --func pbe0
