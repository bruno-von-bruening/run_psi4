from . import *
from run_psi4.types.job import psi4_job

@val_call
def get_single_point(job:psi4_job) -> List[str]:
    
    method=job.method.method
    molecule_name=job.var_names.molecule_name
    storage_var=job.var_names.storage
    jobname=job.jobname

    meth_fl=p4ut.get_method_flags(method)
    wfn_var='wfn'
    if meth_fl['mp2'] or meth_fl['cc']:
        com=['# MP2 and CCSD energies do not correct the density matrix! But the gradient calculation does it!']
        action='gradient'
    else:
        com=[]
        action='gradient'
    if action=='gradient':
        grad_var=f", \'gradient\': {wfn_var}.gradient().to_array().tolist()"
    else:
        grad_var=''
    E_var="E_sp"
    grad_var="Grad"
    lines=com+[
        f"{E_var},{wfn_var} = energy(\'{method}\', molecule={molecule_name}, return_wfn=True)",
        f"{grad_var},{wfn_var} = gradient(\'{method}\', molecule={molecule_name}, return_wfn=True)",
    ]

    if isinstance(jobname, type(None)):
        fchk_file=f"{molecule_name}.fchk"
    else:
        fchk_file=f"{jobname}.fchk"
    lines+=[
        f"# Save the wave function",
        f"fchk_file=\'{fchk_file}\'",
        'fchk(wfn, fchk_file)',
        f"{storage_var}[\'files\'].update({{ \'final_fchk\':fchk_file}})",
        f"{storage_var}[\'results\'].update({{ \'energy\': {E_var}, \'gradient\': {grad_var}.to_array().tolist() }})",
    ]
    return lines
