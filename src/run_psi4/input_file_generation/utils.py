
from . import *
from run_psi4.types.job import psi4_job

ind=4*' '

@val_call
def psi4_action(job,action:str=Field(min_length=1), save_under:str=Field(default='_', min_length=1 ), **kwargs):

    line=save_under
    if 'return_wfn' in kwargs.keys():
        if kwargs['return_wfn']:
            line+=', wfn'

    if not 'molecule' in kwargs.keys():
        kwargs.update(
        molecule=job.var_names.molecule_name,
        )
    

    opts=' , '.join([job.method.func_var]+[f"{k}={v}" for k,v in kwargs.items()] )
    line+=f" = {action}( {opts} )"
    return line


