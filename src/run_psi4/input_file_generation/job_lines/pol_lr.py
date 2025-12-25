from . import *

from run_psi4.types.job import psi4_job, store

from .pol_funcs import center_of_charge_function
from inspect import getsource

def save_fchk(job:psi4_job, wfn_var='wfn'):
    return [
        f"# Save the wave function",
        f"try:",
        ind+f"fchk_file=\'{job.jobname}.fchk\'",
        ind+f"fchk({wfn_var}, fchk_file)",
        ind+f"",
        *[ ind+x for x in  store(job, ('files','final_fchk'), 'fchk_file') ],
        f"except Exception as e:",
        ind+f"print_out(f\'Could not save fchk file due to {{e}}.\')", 
    ]

def prop_calc(job:psi4_job):
    return [
    "# perform polarizability calculation",
    "props = ['DIPOLE_POLARIZABILITIES','QUADRUPOLE_POLARIZABILITIES']",
    f"E, wfn = properties(method, molecule={job.var_names.molecule_name}, properties= props, return_wfn=True)",
    *save_fchk(job, wfn_var='wfn'),
    *store(job, ('results','energetic','energy'), 'E')
    ]

def set_origin(job:psi4_job):
    func_name, func_lines = center_of_charge_function()
    return (func_lines, [
        f"# Define center of charge (in angstrom!) and assert that it has been read in correctly!",
        f"{get_var('center_of_charge')} = {func_name}({job.var_names.molecule_name}) # For evaluation of multipoles and setup of external field",
        f"core.set_global_option('PROPERTIES_ORIGIN', np.array({get_var('center_of_charge')}).tolist() )",
        *store(job, ('results','properties','expansion_center'), f"get_option('None','PROPERTIES_ORIGIN')"),
        # f"assert {get_var('center_of_charge')}==get_option('None','PROPERTIES_ORIGIN'), 'Apparently center of ORIGIN was set incorrectly!'",
    ])
# Get the polarisabilities
def recover_polarisabilities_linrsp(method, the_variables):
    """ Polarizabilities are stored under e.g. 'CCSD(T) DIPOLE POLARIZABILITY XX' etc.
    Need to provide method since it changes the variable polarizabilities are saved under
    Provide the dictionary of psi4
    """
    # Define the keys under which psi4 stores the polarisabilities (important to get cc and mp!)
    key_per_rank={
        'dipole':'DIPOLE POLARIZABILITY',
        'quadrupole':'QUADRUPOLE POLARIZABILITY'
    }
    if any( method.lower().startswith(x) for x in [ 'cc','mp'] ): key_per_rank=dict(
        (k,method.upper()+' '+v) for k,v in key_per_rank.items()
    )

    # Recover all the elements of the polarisability tensor (product of dipole and quadrupole entries)
    dip_entries=['x','y','z']
    quad_entries=['xx','yy','zz','xy','xz','yz']
    keys_to_cover=dict()
    for x in [dip_entries, quad_entries]:
        for k in [
            (a,b) for a in x for b in x
        ]:
            keys_to_cover[k]=False

    # Recover all the entries
    pol_elements=dict()
    for a,b in keys_to_cover.keys():
        full_index=' '.join([key_per_rank['dipole' if len(a)==1 else 'quadrupole'].upper(), a.upper()+b.upper()])
        assert not keys_to_cover[(a,b)], f"Double read in of {{full_index}}?"
        keys_to_cover[(a,b)]=True
        value=the_variables[full_index]
        key=f"{a}@{b}"
        pol_elements.update({key:value})
    assert all(keys_to_cover.values()), f"Not all polarisability entries have been read in: {{keys_to_cover}}"
    return pol_elements


def generate_pol_lr_job_lines(job:psi4_job):
    """ Given a geometry calculate the Polarisabilities available in Psi4
    [ ] Provide starting guess (idk if that will be that much faster though)
    """

    # Get method information

    func_lines=[]

    var_mem.update({'center_of_charge':'center_of_charge'})

    # Polarizabilitiy origin
    func_lines+=set_origin(job)[0] 
    func_lines+=getsource(recover_polarisabilities_linrsp).split('\n')

    pol_store_var='pol_linresp_els'
    lines=[
        f"set_options(dict( e_convergence=1.e-10))",
        # Define the origin (set in psi4 global variables)
        * set_origin(job)[1],
        # Property calculation
        *prop_calc(job),
        # Recover polarisabilities
        f"{pol_store_var}={recover_polarisabilities_linrsp.__name__}({job.method.func_var}, variables())",
        *store(job, ('results','properties','MolPol_LinRsp'), pol_store_var)
    ]

    return lines, func_lines