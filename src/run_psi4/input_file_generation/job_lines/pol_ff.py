from . import *

from run_psi4.types.job import psi4_job, store

from .pol_ff_functions import (
    meassure_multipoles_func, get_moment_function, make_external_field, pol_transform,
    tensor_tag_from_indices, tag_matrix, make_moments
)

lis=List[str]

class options(BaseModel):
    #cc_convergence: float = 1.e-12
    do_moments: bool = True
    delta_dip: float = 1.e-3
    delta_qad: float = 1.e-4
    e_convergence: float|None=None
    def __init__(self, **data):
        super().__init__(**data)
        if self.e_convergence is None:
            self.e_convergence=self.delta_qad**2/100 # Should be 0.01 exact


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
def gather_polarisabilities(job:psi4_job, target_var:str, get_mom=True) -> lis:
    lines=[]

    do_print=True
    if do_print and get_mom:
        lines+= [ f"for t in tensor_elements:",
            ind+f"pol = {get_var('pol_reg')}[t]",
            ind+f"for k in sorted(pol.keys()):",
            ind+ind+f"m=pol[k]",
            ind+ind+f"np.set_printoptions(precision=3, suppress=True)",
            ind+ind+f"the_string=np.array2string(np.array(m),separator = ' ')",
            ind+ind+f"print_out( f\"\\n{{k}}@{{t}}: pol is \\n{{the_string}}\\n\" )",
        ]
    
    # Save
    lines+=[
        f"{target_var}=pol_transform({get_var('pol_reg')})",
    ]
    return lines

@val_call
def generate_pol_ff_job_lines(job:psi4_job, opts:options=options())->Tuple[lis,lis]:
    """Generate job lines to compute polarisability tensors and store them in a yaml file.

    Args:
        job (psi4_job): The psi4 job object.
    """
    lines=[]
    func_lines=[]
    #    ([] if not job.method.method.lower().startswith('cc') else [ ind+f"cc_convergence={opts.cc_convergence},"]) 
    lines+=[
        f"# For HF reference calculation as I understand. Should be probably a bit tighter for finite field (possibly its even faster do to better guess?)",
        f"e_convergence_val={opts.e_convergence}",
        f"set_options(dict(",
        ind+f"e_convergence=e_convergence_val, ",]+[
        f"))",
        *store(job, ('settings','e_convergence'),"e_convergence_val",),
    ]
    # +( [] if not job.method.lower().startswith('cc') else [ f"set_options(dict( cc_convergence=1.e-12 ) )" ]

    func_lines+=make_functions(job)


    # save the external potentials in a dictionary

    from .pol_funcs import center_of_charge_function
    from .pol_ff_functions import make_external_field
    func_lines+=inspect.getsource(make_external_field).split('\n')
    func_name, func_lines_new = center_of_charge_function()
    func_lines+=func_lines_new


    # Get center of charge to define multipoles around
    center_of_charge_var='center_of_charge'
    var_mem.update({'center_of_charge':center_of_charge_var})
    lines+=[ 
        f"# External FIELD PARAMETERS",
        f"{get_var('center_of_charge')} = {func_name}({job.var_names.molecule_name}) # For evaluation of multipoles and setup of external field",
        f"distance=1000.0 # Large distance to approximate uniform field" ,
        f"delta_dip={opts.delta_dip:E}",
        f"delta_qad={opts.delta_qad:E}",
        *store(job, ('settings','ff_val_dip'),'delta_dip'),
        *store(job, ('settings','ff_val_qad'),'delta_qad'),
    ]
    if len(job.geom.atom_types)==1:
        com=" # Single atom only need these elements"
        tags=['x','xx','xy']
    else:
        tags=['x','y','z','xx','xy','yy','xz','yz','zz']#,'xy','yz','xz']
        com=''
    lines+=[ 
        f"",
        f"# Field moments to apply finited field",
        f"tensor_elements=[{com}" ]+[
            ind+f"'{tag}'," for tag in tags ]+[
        f"]"
    ]

    lines+=[
        f"",
        f"psi4.set_global_option('PROPERTIES_ORIGIN', np.array({get_var('center_of_charge')}).tolist() )",
    ]
    lines+=[ 
        f"# vec=[ -float(x) for x in {get_var('center_of_charge')} ]",
        f"# {job.var_names.molecule_name}.translate(core.Vector3(vec[0],vec[1],vec[2]))",
    ] 
    # lines+=[f"psi4.set_options(dict(  MOMENTS_ORIGIN={get_var('center_of_charge')}.tolist()  ))"
    #     ]

    var_mem.update({'pot_registry':'pot_reg'}) 
    lines+=[
        f"",
        f"# Generate external potentials (unit value will be tuned by FF value)",
        f"{get_var('pot_registry')}=dict() # Save all external potentials there (base version)",
        f"for tensor_key in tensor_elements:",
        ind+f"pot_base = make_external_field(tensor_key, {get_var('center_of_charge')}, d=distance)",
        ind+f"{get_var('pot_registry')}.update({{ tensor_key:pot_base }})",
        f"{get_var('pot_registry')}.update( zero=None ) # zero potenial",
    ]

    # Evaluate moments for each
    do_single_point_mom=job.method.can_do_gradient
    var_mem.update({f"ff_reg": 'ff_reg'}),
    lines+=[
        f"",
        f"# Evaluate wfn in external potential",
        f"",
        f"{get_var('ff_reg')}=dict()",
        f"E,wfn,grad,mom= get_moments( {get_var('pot_registry')}['zero'], 0, do_mom={do_single_point_mom} ) # Should be True if we want zero point moments",
        f"{get_var('ff_reg')}.update( zero=dict({{ 'delta':0, 'E':E, 'mom': mom }}))",
        f"",
        f"# Save the wave function",
        f"fchk_file=\'{job.jobname}.fchk\'",
        'fchk(wfn, fchk_file)',
        f"",
        f"if not os.path.isfile(fchk_file):",
        ind+f"fchk_file=None",
        f"{job.var_names.storage}[\'files\'].update({{ \'final_fchk\':fchk_file}})",
        f"{job.var_names.storage}[\'results\']['energetic'].update({{ \'energy\': E, \'gradient\': (grad.to_array().tolist() if grad is not None else None) }})",
        f"for i,t in enumerate(tensor_elements):",
        ind+f"delta= (delta_dip if len(t)==1 else delta_qad)",
        ind+f"for x in [ -1 , 1 ]: # negative and positive field",]+ \
        [ind+ind+f"do_mom=True if x>0 else False  # only do moments for positive field (negative is for energy derviation) "
          if opts.do_moments else ind+ind+f"do_mom=False  #no moments request by main script"
        ]+[
        ind+ind+f"E,wfn,grad,mom= get_moments( {get_var('pot_registry')}[t] , x*delta, ",]+\
        ([3*ind+f"do_mom=False ) #no moments request by main script"] if not opts.do_moments else [
        ind+ind+ind+f"do_mom=do_mom ) # only do moments for positive field (negative is for energy derviation) ",])+[
        ind+ind+f"{get_var('ff_reg')}.update({{ t+('_neg' if x==-1 else '_pos'):dict({{ 'delta': x*delta, 'E': E, 'mom': mom }} )}})",
        ind+f"print_out(f\"CUSTOM MARKER BY USER: Finished moment '{{t}}' ({{i}} out of {{len(tensor_elements)}}, {{i+1/len(tensor_elements)*100}}%)\")",
    ]


    # Compute the finite difference
    var_mem.update({f"pol_reg":"pol_reg"})
    lines += [ 
        f"",
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
        ind+ind+f"*( 1/2 if len(t)==2 else 1) ) # Induction Energy",
        ind+f"if any(x['mom'] is None for x in [first, second] ): continue",
        ind+f"for k in set(list(first['mom'].keys())+list(second['mom'].keys()) ):",
        ind+ind+f"assert k in first['mom'].keys() and k in second['mom'].keys()",
        ind+ind+f"ff_mom=( first['mom'][k] - second['mom'][k]) / delta",
        ind+ind+f"cur_pol.update({{ k:  np.array(ff_mom).tolist() }})",
        ind+f"{get_var('pol_reg')}.update({{ t : cur_pol }})",
    ]


    lines+=[
        f"",
        f"# Gather and reform polarisabilities",
    ]+gather_polarisabilities(job, "pol_vals", get_mom=opts.do_moments)
    


    path=('results','properties')
    lines+=[f"#Store object",  
            *store(job,path+('expansion_center',), f"[ float(x) for x in {get_var('center_of_charge')}]"),
            *store(job,path+('MolPol_FinFie_through_dens',), 'pol_vals'),
            *store(job,path+('MolPol_FinFie_through_eng',), 'pol_through_energy') ,
    ]
    if do_single_point_mom:
        lines+=make_moments(f"{get_var('ff_reg')}['zero']['mom']",'moms')
        lines+=store(job,path+('MolMom',), 'moms')


    return lines, func_lines



# set_options( {'perturb_dipole': make_c2v_vec(delta),'reference':'rhf' })
# set {
# OPDM_RELAX True
# }
