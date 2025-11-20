#!/usr/bin/env python

import os, sys, re, time
import numpy as np

from run_psi4.main.main_base import psi4_run

def main():
    description=None
    epilog=None
    job_choices=['polarizabilities','single_point']
    import argparse as ap; par=ap.ArgumentParser(description=description, epilog=epilog)
    add=par.add_argument
    add('job', help=f"What to do", choices=job_choices)
    add('--xyz', help=f"xyz_file to process")
    add('--func', help=f"Functional for psi4 calculations")
    add('--basis', help=f"Main Basis set for Psi4 Calculation")
    add('--exc', action='store_true', help=f"Exceute psi4 (requires conda environment that runs psi4)")
    add('--test', action='store_true', help=f"Run a test (also provide exc statement to run!)")
    add(
        '--num_threads', type=int, help=f"Number of threads the psi4 process will use", default=8,
    )
    add(
        '--memory', type=str, help=f"Memory size (has to end with MB or GB)", default='8GB',
    )
    args=par.parse_args()
    xyz_file    = args.xyz
    job         = args.job
    func        = args.func
    basis_set   = args.basis
    execute     = args.exc
    do_test     = args.test
    num_threads = args.num_threads
    memory      = args.memory

        

    psi4_run(
        job,
        func=func,
        basis_set=basis_set,
        xyz_file=xyz_file,
        execute=execute,
        do_test=do_test,
        num_threads=num_threads,
        memory=memory,
    )

if __name__=='__main__':
    main()

