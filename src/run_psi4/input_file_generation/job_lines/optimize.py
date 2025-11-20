from . import *
from run_psi4.types.job import psi4_job

@val_call
def optimization(job:psi4_job) -> List[str]:
    """ We request the same as for single point but additionally with optimization instead of energy"""
    
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

    kwargs=dict(
        return_wfn=True,
        molecule=molecule_name
    )
    lines=com+[
        psi4_action(job, 'optimize', save_under=E_var, **kwargs),
        psi4_action(job, 'gradient', save_under=grad_var, **kwargs),
    ]

    if isinstance(jobname, type(None)):
        stem=f"{molecule_name}"
    else:
        stem=f"{jobname}"
    job.files.fchk=stem+".fchk"

    lines+=[
        f"# Save the wave function",
        f"fchk_file=\'{job.files.fchk}\'",
        'fchk(wfn, fchk_file)',
        f"try:",
        ind+f"xyz_geom=Molecule.from_string({job.var_names.molecule_name}.save_string_xyz() ).save_string_xyz()",
        #ind+f"{job.var_names.molecule_name}.reorient()",
        f"except Exception as ex:",
        ind+f"xyz_geom={job.var_names.molecule_name}.save_string_xyz()",
        ind+f"print_out(f\"WARNING: could not cannonically orient the molecule: {{str(ex)}}\")",
    ]

    saveargs=dict(
        files=dict(
            final_fchk='fchk_file',
        ),
        results=dict(
            energy=E_var,
            gradient=f"np.array({grad_var}, dtype=float).tolist()",
            optimized_geometry=f"xyz_geom",
        )
    )

    lines+=update_storage(job,saveargs)

    return lines