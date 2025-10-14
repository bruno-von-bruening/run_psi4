#!/usr/bin/env python

import os, sys, re, time
import numpy as np

from run_psi4.main.main import psi4_run

def main():
    description=None
    epilog=None
    job_choices=['polarizabilities','single_point']
    import argparse as ap; par=ap.ArgumentParser(description=description, epilog=epilog)
    add=par.add_argument
    
    add('job', help=f"What to do", choices=job_choices)
    
    # System specification
    add('--xyz', help=f"xyz_file to process")
    add('--total_charge', type=int,help=f"Charge of the system", default=0)
    add('--multiplicity', type=int, default=1, help=f"Multiplicity of system spin state")
    add('--ac_shift', type=float, default=None, help=f"AC shift")
    
    # Calculation specification
    add('--func', default='pbe0', help=f"Functional for psi4 calculations")
    add('--basis', default='aug-cc-pVDZ', help=f"Main Basis set for Psi4 Calculation")
    
    # Run specifications
    add('--exc', action='store_true', help=f"Exceute psi4 (requires conda environment that runs psi4)")
    add('--executable', help=f"Exceutable to run psi4 with!")
    add('--test', action='store_true', help=f"Run a test (also provide exc statement to run!)")
    add(
        '--num_threads', type=int, help=f"Number of threads the psi4 process will use", default=8,
    )
    add(
        '--memory', type=str, help=f"Memory size (has to end with MB or GB)", default='8GB',
    )
    args=par.parse_args()
    xyz_file    = args.xyz
    total_charge= args.total_charge
    multiplicity= args.multiplicity
    ac_shift    = args.ac_shift
    # 
    job         = args.job
    func        = args.func
    basis_set   = args.basis
    execute     = args.exc
    executable  = args.executable
    do_test     = args.test
    num_threads = args.num_threads
    memory      = args.memory

    psi4_run(
        job,
        func=func,
        basis_set=basis_set,
        total_charge=total_charge,
        multiplicity=multiplicity,
        xyz_file=xyz_file,
        execute=execute,
        executable=executable,
        do_test=do_test,
        num_threads=num_threads,
        memory=memory,
        ac_shift=ac_shift,
    )
                

if __name__=='__main__':
    main()
