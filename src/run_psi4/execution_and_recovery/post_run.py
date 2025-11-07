from . import *

class exit_status(Enum):
    normal='normal'
    failed='failed'

@val_call
def psi4_exit_ana(output_file:pdtc_file) -> exit_status:
    """ CHECK that all went well """

    with open(output_file, 'r') as rd:
        lines=rd.readlines()
    last_line=lines[-1]
    successful_tag=r'\*\*\* Psi4 exiting successfully[.]*^'
    
    normal_exit_line_fl=bool(
        re.match(successful_tag, last_line)
    )

    if normal_exit_line_fl:
        exit_status=exit_status.normal
    else:
        exit_status=exit_status.failed

    return exit_status