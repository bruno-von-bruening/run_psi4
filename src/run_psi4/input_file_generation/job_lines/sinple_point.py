from . import *
from run_psi4.types.job import psi4_job
from run_psi4.types.method import method_types


from run_psi4.input_file_generation.functions.get_moments import handle_molecular_moments

@val_call
def get_single_point(job:psi4_job) -> List[str]:
    
    method=job.method.method
    molecule_name=job.var_names.molecule_name
    storage_var=job.var_names.storage
    jobname=job.jobname

    meth_fl=p4ut.get_method_flags(method)
    wfn_var='wfn'

    han_mom=handle_molecular_moments()
    set_orig_func,orig_lines=han_mom.set_origin(job)
    

    func_lines=han_mom.make_functions(job)+set_orig_func

    if not job.method.can_do_gradient:
        raise ValueError(f"Cannot do gradient for {job.method}")

    if job.method.method_type in [ method_types.MP, method_types.CC ]:
        com=['# MP2 and CCSD energies do not correct the density matrix! But the gradient calculation does it!']
        action='gradient'
    else:
        com=[]
        action='gradient'
    if action=='gradient':
        grad_var=f", \'gradient\': {wfn_var}.gradient().to_array().tolist()"
    else:
        grad_var=''

    # raise Exception(action, job.method)
    E_var="E_sp"
    grad_var="Grad"
    lines=com+func_lines+[
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

    mol_mom_var=f"mol_moms"

    lines+=[ 
        *orig_lines,
        "",
        *han_mom.eval_moments_to(job, mol_mom_var), 
        f"{storage_var}[\'results\'].update({{ \'MolMoms_from_dens\': {mol_mom_var}}})" ]
    return lines
