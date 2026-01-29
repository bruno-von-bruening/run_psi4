#!/usr/bin/env python

# Need to have this package installed! (side project from qcp_global_utils project)
from qcp_versioning.package_docu import main_parser



MAIN_PKG='run_psi4'
DEPENDENCIES=[
    'qcp_global_utils',
    'qcp_objects',
]

if __name__ == '__main__':
    main_parser(MAIN_PKG, DEPENDENCIES)