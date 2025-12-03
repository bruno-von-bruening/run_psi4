from . import *

from run_psi4.types.job import psi4_job

from .pol_ff_functions import (
    meassure_multipoles_func, get_moment_function, make_external_field, pol_transform,
    tensor_tag_from_indices, tag_matrix
)

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
    for x in [tensor_tag_from_indices, tag_matrix, pol_transform ]:
        lines+=f"{inspect.getsource(x)}".split('\n')

    return lines

@val_call
def gather_polarisabilities(job:psi4_job) -> lis:
    lines=[]

    do_print=True
    if do_print:
        lines+= [ f"for t in tensor_elements:",
            ind+f"pol = {get_var('pol_reg')}[t]",
            ind+f"for k in sorted(pol.keys()):",
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
    lines+=[
        f"set_options(dict( e_convergence=1.e-12 )) # For HF reference calculation as I understand. Should be probably a bit tighter for finite field (possibly its even faster do to better guess?)"
    ]
    # +( [] if not job.method.lower().startswith('cc') else [ f"set_options(dict( cc_convergence=1.e-12 ) )" ]

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

    from .pol_ff_functions import center_of_charge_function


    # Get center of charge to define multipoles around
    func_name, new_lines = center_of_charge_function()
    lines+=new_lines
    center_of_charge_var='center_of_charge'
    var_mem.update({'center_of_charge':center_of_charge_var})
    lines+=[
        f"{get_var('center_of_charge')} = {func_name}({job.var_names.molecule_name})",
    ]

    from .pol_ff_functions import make_external_field
    lines+=inspect.getsource(make_external_field).split('\n')
    tags=['x','y','z','xx','yy','zz']#,'xy','yz','xz']
    lines+=[ 
        f"tensor_elements=[" ]+[
            ind+f"'{tag}'," for tag in tags ]+[
        f"]"
    ]
    var_mem.update({'pot_registry':'pot_reg'}) 
    lines+=[ f"distance=1000.0 # Large distance to approximate uniform field" ]
    lines+=[f"{get_var('pot_registry')}=dict() # Save all external potentials there (base version)"]
    lines+=[
        f"for tensor_key in tensor_elements:",
        ind+f"pot_base = make_external_field(tensor_key, {get_var('center_of_charge')}, d=distance)",
        ind+f"{get_var('pot_registry')}.update({{ tensor_key:pot_base }})"
    ]



    # Evaluate moments for each
    var_mem.update({f"ff_reg": 'ff_reg'}),
    lines+=[
        f"",
        f"# Tune potentials (give finite difference value)",
        f"{get_var('pot_registry')}.update( zero=None )", # Define zero potential
        f"delta_dip={1.e-4:E}",
        f"delta_qad={1.e-5:E}",
        f"",
        f"{get_var('ff_reg')}=dict()",
        f"E,wfn,grad,mom= get_moments( {get_var('pot_registry')}['zero'], 0 )",
        f"{get_var('ff_reg')}.update( zero=dict({{ 'delta':0, 'E':E, 'mom': mom }}))",
        f"# Save the wave function",
        f"fchk_file=\'{job.jobname}.fchk\'",
        'fchk(wfn, fchk_file)',
        f"{job.var_names.storage}[\'files\'].update({{ \'final_fchk\':fchk_file}})",
        f"{job.var_names.storage}[\'results\'].update({{ \'energy\': E, \'gradient\': (grad.to_array().tolist() if grad is not None else None) }})",
        f"for t in tensor_elements:",
        ind+f"delta= (delta_dip if len(t)==1 else delta_qad)",
        ind+f"for x in [ -1 , 1 ]: # negative and positive field",
        ind+ind+f"E,wfn,grad,mom= get_moments( {get_var('pot_registry')}[t] , x*delta, ",
        ind+ind+ind+f"do_mom=True if x>0 else False ) # only do moments for positive field (negative is for energy derviation) ",
        ind+ind+f"{get_var('ff_reg')}.update({{ t+('_neg' if x==-1 else '_pos'):dict({{ 'delta': x*delta, 'E': E, 'mom': mom }} )}})"
    ]


    # Compute the finite difference
    var_mem.update({f"pol_reg":"pol_reg"})
    lines += [ 
        f"# get polarisabilities through finite differences",
        f"{get_var('pol_reg')}=dict()",
        f"pol_through_energy=dict()",
        f"for t in tensor_elements:",
        ind+f"first,second, third= {get_var('ff_reg')}[t+'_pos'], {get_var('ff_reg')}['zero'], {get_var('ff_reg')}[t+'_neg']",
        ind+f"print_out(f\"E={{first['E']}} , E0={{second['E']}} for tag={{t}} and delta={{delta}}\\n\")",
        ind+f"cur_pol=dict()",
        ind+f"delta = (first['delta'] - second['delta'] )",
        ind+f"mapper=dict( (k,i) for i,k in enumerate(['x','y','z']) )",
        ind+f"pol_through_energy[f\"{{t}}@{{t}}\"]=float(-(first['E']-2*second['E']+third['E'])/(delta**2) \ # FD formula (polarisabilities are positive by convention)",
        ind+ind+f"*( 1/2 if len(t)==2 else 1) ) # Do not know where this comes from (convention in psi4?)",
        ind+f"for k in set(list(first['mom'].keys())+list(second['mom'].keys()) ):",
        ind+ind+f"assert k in first['mom'].keys() and k in second['mom'].keys()",
        ind+ind+f"ff_mom=( first['mom'][k] - second['mom'][k]) / delta",
        ind+ind+f"cur_pol.update({{ k:  np.array(ff_mom).tolist() }})",
        ind+f"{get_var('pol_reg')}.update({{ t : cur_pol }})",
    ]

    lines+=[
        f"",
        f"# Gather and reform polarisabilities",
    ]+gather_polarisabilities(job)+[
        f"{job.var_names.storage}['results'].update( pol_ff_through_energy = pol_through_energy )",
    ]

    lines+=[
        f"# Store molecular Multipoles",
        f"moms={{}}",
        f"for k,v in {get_var('ff_reg')}['zero']['mom'].items():",
        f"    moms.update({str(tag_matrix.__name__)}(v))",
        f"{job.var_names.storage}['results'].update( multipoles = moms )",
    ]

    return lines



# set_options( {'perturb_dipole': make_c2v_vec(delta),'reference':'rhf' })
# set {
# OPDM_RELAX True
# }
