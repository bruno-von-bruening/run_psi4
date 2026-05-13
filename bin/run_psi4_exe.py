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
    epilog='\n'.join([
        f"Usage examples:",*[ f"   -  {k} : {v}" for k,v in (
        (f"Get a example config file for given job", f"{os.path.basename(__file__)} --job_type <select_job_type>"),
        )]
    ])
    import argparse as ap; par=ap.ArgumentParser(description=description, epilog=epilog, formatter_class=ap.RawDescriptionHelpFormatter)
    add=par.add_argument


    add(  'config',     
        nargs='?', type=str, default=None, help=f"Config in yaml or json format, default one (that will be used if executed) will be dropped if non is provided")
    add(  '--exc', '-x',
        action='store_true', help=f"Exceute psi4 (requires conda environment that runs psi4)")
    add(  '--executable',
        type=str, help=f"Exceutable to run psi4 with!")
    add(  '--job_type', '-j', 
        help=f"What to do", choices=allowed_job_choices)
    add(  '--method', '-m', 
        help=f"Method to use (overrides config file)", default=None)
    add(  '--basis', '-b', 
        help=f"Basis to use (overrides config file)", default=None)
    add(  '--verbose', '-v', 
        action='store_true', help=f"Verbose output")
    add(  '--edit', '-e',
         action='store_true', help=f"Open the config file in vim after creating it (only when creating example config file)")
    add(  '--geom', '-g', 
        type=str, help=f"Geometry to include in input file")
    add( '--version_info','-vinfo', action='store_true')
    args=par.parse_args()

    print_version_info=args.version_info
    if print_version_info:
        from run_psi4._version import __commit_id__
        import run_psi4
        print(f"{run_psi4.__name__} {__commit_id__}")
        return None


    if args.verbose:
        logger.setLevel(VERBOSE_LEVEL_NUM)
    else:
        logger.setLevel(logging.INFO)  # default; set to VERBOSE_LEVEL_NUM or logging.DEBUG for verbose output

    config_file=args.config
    exc=args.exc
    executable=args.executable
    job_type=args.job_type
    edit=args.edit

    method=args.method
    geom=args.geom
    basis=args.basis
    #method=args.method

    if executable:
        if not os.path.isfile(executable): raise FileNotFoundError(f"Provided exectuable \'{executable}\' is not a file!")

    make_example= config_file is None

    if not make_example:
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
    else:
        config_file=_config_file_default
        if not job_type:
            # Show usage and quit
            par.print_help()
            quit(f"\n\nError: Did not provide job_type through command line or config file! (see help above)")
        kwargs={}
        if any( [ x is not None for x in (method, basis)]):
            kwargs.update( method=dict(
                method_tag=method,
                basis_set=basis
            ))
        if geom:
            if geom.endswith('.xyz'):
                kwargs['geom']=dict(
                    xyz_file=geom,
                    charge=0,
                    multiplicity=1,
                )

        config=psi4_job.example(execute=exc, executable=executable, jobtag=job_type, **kwargs)

        with open(config_file,'w') as wr:
            yaml.safe_dump(config, wr)

    if edit:
        try:
            os.system(f"nvim {config_file}")
        except Exception as e:
            warn(f"Could not open {config_file} in nvim: {e}")
        config=load_json_or_yaml(config_file)
    
    try:
        job=psi4_job(**config)
    except Exception as e:
        raise RuntimeError(f"Config file {config_file} does not translate into valid psi4_job: {e}")

    psi4_run(job)

if __name__=='__main__':
    main()
