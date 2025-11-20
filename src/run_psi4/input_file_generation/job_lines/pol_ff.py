from . import *

from run_psi4.types.job import psi4_job

from .pol_ff_functions import meassure_multipoles_func, get_moment_function, make_external_field, pol_transform

lis=List[str]

@val_call
def make_functions(job:psi4_job):
    # make functions, that will be used later
    lines=[
        *meassure_multipoles_func(job),
        '',
        *get_moment_function(job),
        '',
    ]
    lines+=f"{inspect.getsource(pol_transform)}".split('\n')

    return lines

@val_call
def gather_polarisabilities(job:psi4_job) -> lis:
    lines=[]

    do_print=True
    if do_print:
        lines+= [ f"for t in tags:",
            ind+f"pol = {get_var('pol_reg')}[t]",
            ind+f"for k in pol.keys():",
            ind+ind+f"m=pol[k]",
            ind+ind+f"np.set_printoptions(precision=3, suppress=True)",
            ind+ind+f"the_string=np.array2string(np.array(m),separator = ' ')",
            ind+ind+f"print_out( f\"\\n{{k}}@{{t}}: pol is \\n{{the_string}}\\n\" )",
        ]
    
    # Save
    storage_var=job.var_names.storage
    lines+=[
        f"pol_trans=pol_transform({get_var('pol_reg')})",
        f"{storage_var}.update( pol_ff = pol_trans)"
    ]
    return lines

@val_call
def generate_pol_ff_job_lines(job:psi4_job)->lis:
    """Generate job lines to compute polarisability tensors and store them in a yaml file.

    Args:
        job (psi4_job): The psi4 job object.
    """
    lines=[]

    function_lines=make_functions(job)
    lines+=function_lines
    # file=f"functions_for_pol_ff.py"
    # with open(file, 'w') as wr:
    #     wr.write('\n'.join(function_lines))
    
    # lines+=[ 
    #     f"sys.path.insert(1,__file__)",
    #     f"raise Exception(os.path.dirname(__file__))",
    #     f"from {file.replace('.py','')} import *"
    # ]


    # save the external potentials in a dictionary
    var_mem.update({'pot_registry':'pot_reg'}) 
    lines+=[f"{get_var('pot_registry')}=dict() # Save all external potentials there (base version)"]

    tags=['x','y','z','xx','yy','zz','xy','yz','xz']
    for tag in tags:
        lines +=make_external_field(tag) # makes external corresponding to perturbation and adds to registry

    # Evaluate moments for each
    var_mem.update({f"ff_reg": 'ff_reg'}),
    lines+=[
        f"",
        f"# Tune potentials (give finite difference value)",
        f"{get_var('pot_registry')}.update( zero=None )", # Define zero potential
        f"tags={tags}",
        f"delta={1.e-8}",
        f"",
        f"{get_var('ff_reg')}=dict()",
        f"{get_var('ff_reg')}.update( zero=dict({{ 'delta':0, 'mom': get_moments({get_var('pot_registry')}['zero'], 0)}} ))",
        f"for t in tags:",
        ind+f"{get_var('ff_reg')}.update({{ t:dict({{ 'delta': delta, 'mom': get_moments({get_var('pot_registry')}[t] , delta) }} )}})"
    ]  


    # Compute the finite difference
    var_mem.update({f"pol_reg":"pol_reg"})
    lines += [ 
        f"# get polarisabilities through finite differences",
        f"{get_var('pol_reg')}=dict()",
        f"for t in tags:",
        ind+f"first,second= {get_var('ff_reg')}[t], {get_var('ff_reg')}['zero']",
        ind+f"cur_pol=dict()",
        ind+f"for k in set(list(first['mom'].keys())+list(second['mom'].keys()) ):",
        ind+ind+f"assert k in first['mom'].keys() and k in second['mom'].keys()",
        ind+ind+f"ff_mom=( first['mom'][k] - second['mom'][k]) / (first['delta'] - second['delta'] )",
        ind+ind+f"cur_pol.update({{ k:  np.array(ff_mom).tolist() }})",
        ind+f"{get_var('pol_reg')}.update({{ t : cur_pol }})",
    ]

    lines+=[
        f"",
        f"# Gather and reform polarisabilities",
    ]+gather_polarisabilities(job)

    return lines



# set_options( {'perturb_dipole': make_c2v_vec(delta),'reference':'rhf' })
# set {
# OPDM_RELAX True
# }
