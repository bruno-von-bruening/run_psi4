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
    if any( job.method.method.lower().startswith(x) for x in ['cc','mp']):
        com=['# MP2 and CCSD energies do not correct the density matrix! But the gradient calculation does it!']
    else: com=[]
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
    opt_kwargs=dict(
        geom_maxiter=job.job_settings.maxiter,
        # g_convergence=job.job_settings.g_convergence,
        # e_convergence=job.job_settings.e_convergence,
    )

    do_gradient=True

    lines=com+[
        f"# Settings specific to geometry optimization",
        f"psi4.set_options(dict(", 
        '\n'.join([f"    {key}={value}," for key, value in opt_kwargs.items()]), 
        f"))",
        f"",f"# Perform the geometry optimization",
        f"start_opt=time.time()",
        *psi4_action(job, 'optimize', save_under=E_var, **kwargs),
        *store_timing(job, 'geometry_optimization', 'start_opt',add=True),
    ]
    if do_gradient:
        lines+=[
            f"",f"# Compute the gradient at the optimized geometry",
            f"start_grad=time.time()",
            *psi4_action(job, 'gradient', save_under=grad_var, **kwargs),
            *store_timing(job, 'comp_grad_after_opt', 'start_grad',add=True),
        ]
    def check_hessian_lines():
        """ Confirm that the Hessian has only positive eigenvalues and the found geometry is indeed a minimum
        Probably there is more to do with that (save it)"""
        # H_var="Hess"
        lin=[   
            *psi4_action(job, 'freq', save_under="E", **kwargs),
            # f"{H_var} = wfn.frequencies()",
            # f"eigvals = np.linalg.eigvalsh({H_var})",
            f"eigvals = np.array(wfn.frequencies().to_array())",
            f"if min(eigvals) <= 1.e-4: raise Exception(f\"Optimization did not reach a minimum, found not positive Hessian eigenvalue: {{min(eigvals)}}\")",
            f"else: print_out(f\"MAKRER: Hessian check passed, all eigenvalues positive, min eigenvalue: {{min(eigvals)}}\")",
        ]
        return lin
    if job.job_settings.check_hessian:                                                          
        lines+= [
            f"",
            f"# Perform a frequency calculation to check the Hessian at the optimized geometry",
            f"start_freq=time.time()",
            *check_hessian_lines(),
            *store_timing(job, 'frequency_check', 'start_freq',add=True)
        ]

    if isinstance(jobname, type(None)):
        stem=f"{molecule_name}"
    else:
        stem=f"{jobname}"
    job.files.fchk=stem+".fchk"

    lines+=[
        f"",
        f"# Save the wave function",
        f"fchk_file=\'{job.files.fchk}\'",
        'fchk(wfn, fchk_file)',
        f"try:",
        ind+f"xyz_geom=Molecule.from_string({job.var_names.molecule_name}.save_string_xyz() ).save_string_xyz()",
        #ind+f"{job.var_names.molecule_name}.reorient()",
        f"except Exception as ex:",
        ind+f"xyz_geom={job.var_names.molecule_name}.save_string_xyz()",
        ind+f"print_out(f\"WARNING: could not cannonically orient the molecule: {{str(ex)}}\")",
        f"with open(\'{job.jobname}.xyz\', \'w\') as wr:",
        ind+f"wr.write(f\"{{{job.var_names.molecule_name}.natom()}}\\nOptimized with psi4;  Char Mult{{xyz_geom}}\")",
    ]

    saveargs=dict(
        files=dict(
            final_fchk='fchk_file',
        ),
        results=dict(
            energy=E_var,
            optimized_geometry=f"xyz_geom",
        )
    )
    if do_gradient:
        saveargs['results']['gradient'] = f"np.array({grad_var}, dtype=float).tolist()"

    lines+=update_storage(job,saveargs)

    return lines
