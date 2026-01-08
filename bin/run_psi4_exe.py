#!/usr/bin/env python

# The name of this file needs to be different that the name of the module.... (would prefer to just call it through run_psi4.py)

import os, sys, re, time
import numpy as np, yaml

from run_psi4.main.main import psi4_run, allowed_job_choices
from run_psi4.util.import_helper import warn, info
from qcp_global_utils.environment.file_handling import load_json_or_yaml

from run_psi4.types.job import psi4_job

GEOM_KEY='geom'
_config_file_default='run_psi4_config.yaml'


from run_psi4.util.import_helper import logger, logging, VERBOSE_LEVEL_NUM

def main():
    description=None
    epilog=None
    import argparse as ap; par=ap.ArgumentParser(description=description, epilog=epilog)
    add=par.add_argument


    add(  'config', nargs='?', type=str, default=None, help=f"Config in yaml or json format, default one (that will be used if executed) will be dropped if non is provided")
    add(  '--exc', action='store_true', help=f"Exceute psi4 (requires conda environment that runs psi4)")
    add(  '--executable', type=str, help=f"Exceutable to run psi4 with!")
    add(  '--job_type', '-j', help=f"What to do", choices=allowed_job_choices)
    add(  '--method', '-m', help=f"Method to use (overrides config file)", default=None)
    add(  '--verbose', '-v', action='store_true', help=f"Verbose output")

    args=par.parse_args()

    if args.verbose:
        logger.setLevel(VERBOSE_LEVEL_NUM)
    else:
        logger.setLevel(logging.INFO)  # default; set to VERBOSE_LEVEL_NUM or logging.DEBUG for verbose output

    config_file=args.config
    exc=args.exc
    executable=args.executable
    job_type=args.job_type
    #method=args.method

    if executable:
        if not os.path.isfile(executable): raise FileNotFoundError(f"Provided exectuable \'{executable}\' is not a file!")

    if config_file:
        assert os.path.isfile(config_file), f"Config file {config_file} not found!"
        config=load_json_or_yaml(config_file)
        if not isinstance(config, dict): raise ValueError( f"Config file {config_file} did not load to a dictionary!" )

        if 'execute' in config:
            warn(f"Overriding execute flag in config file {config_file} with command line argument!")
        config['execute']=exc

        if executable:
            if 'executable' in config:
                warn(f"Overriding executable in config file {config_file} with command line argument!")
            config['executable']=executable

        if 'jobtag' in config:
            if job_type:
                warn(f"Overriding jobtag in config file {config_file} with command line argument!")
                config['jobtag']=job_type
        else:
            if job_type is None:
                raise ValueError(f"Did not provide job_type through command line or config file!")

        if args.method:
            if 'method' in config:
                warn(f"Overriding method in config file {config_file} with command line argument!")
            config['method']=args.method


        # try:
        #     psi4_job(**config)
        # except Exception as e:
        #     raise RuntimeError(f"Config file {config_file} is not valid for psi4_job: {e}")
        
    else:
        config_file=_config_file_default
        if not job_type:
            raise ValueError(f"Did not provide job_type through command line or config file!")
        config=psi4_job.example(execute=exc, executable=executable, jobtag=job_type)

        with open(config_file,'w') as wr:
            yaml.safe_dump(config, wr)
    
    try:
        job=psi4_job(**config)
        
    except Exception as e:
        raise RuntimeError(f"Config file {config_file} does not translate into valid psi4_job: {e}")
    psi4_run(job)

    return


    if False:

        # System specification
        add(f"--{GEOM_KEY}", help=f"geometry to process (currently only xyz file supported)")
        add('--total_charge', type=int,help=f"Charge of the system", default=0)
        add('--multiplicity', type=int, default=1, help=f"Multiplicity of system spin state")
        add('--ac_shift', type=float, default=None, help=f"AC shift")
        
        # Calculation specification
        add('--method', default='pbe0', help=f"Functional for psi4 calculations")
        add('--basis', default='aug-cc-pVDZ', help=f"Main Basis set for Psi4 Calculation")
        
        # Run specifications
        add('--test', action='store_true', help=f"Run a test (also provide exc statement to run!)")
        add(
            '--num_threads', type=int, help=f"Number of threads the psi4 process will use", default=8,
        )
        add(
            '--memory', type=str, help=f"Memory size (has to end with MB or GB)", default='8GB',
        )

        add(
            '--no_mom_ff', action='store_true', help=f"Whether to directly compute and store multipole moments during finite field calculations",
        )
        add(
            '--no_freeze_core','-ncf',action='store_true', help=f"Do not freeze core orbitals in correlated calculations",
        )
        add(
            '--no_df','-ndf',action='store_true', help=f"Do not use density fitting in correlated calculations",
        )
        add(
            '--unrestricted','-ur',action='store_true', help=f"Use unrestricted orbitals for closed shell systems",
        )

        args=par.parse_args()
        xyz_file    = getattr(args, GEOM_KEY)
        total_charge= args.total_charge
        multiplicity= args.multiplicity
        ac_shift    = args.ac_shift
        # 
        job         = args.job
        func        = args.method
        basis_set   = args.basis
        execute     = args.exc
        executable  = args.executable
        do_test     = args.test
        num_threads = args.num_threads
        memory      = args.memory
        do_mom=not args.no_mom_ff
        freeze_core=not args.no_freeze_core
        use_df=not args.no_df
        unrestricted=args.unrestricted

        if total_charge is None:
            warn(f"No total charge provided, assuming neutral system!")
            total_charge=0
        if multiplicity is None:
            warn(f"No multiplicity provided, assuming singlet!")
            multiplicity=1

        if do_test and xyz_file is None:
            geom=None
        else:
            geom= {
                'input':xyz_file,
                'charge': total_charge,
                'multiplicity':multiplicity,
                'ac_shift': ac_shift,
            }

    psi4_run(
        job,
        method=func,
        basis_set=basis_set,
        geom=geom,
        execute=execute,
        executable=executable,
        do_test=do_test,
        machine_settings=dict(
            num_threads=num_threads,
            memory=memory,
        ),
        do_mom=do_mom,
        freeze_core=freeze_core,
        use_df=use_df,
        unrestricted=unrestricted
    )

if __name__=='__main__':
    main()
