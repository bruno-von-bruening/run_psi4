from . import *

from run_psi4.types.job import psi4_job, store

from .pol_ff_functions import (
    meassure_multipoles_func, get_moment_function, make_external_field, pol_transform,
    tensor_tag_from_indices, tag_matrix, make_moments
)

lis=List[str]

from run_psi4.types.job_settings import pol_ff_options



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

def get_second_derivative(plus, minus, zero,delta):
    return (plus-2*zero+minus)/delta**2 
def get_first_derivative(plus, minus, delta):
    return (plus-minus)/(delta) 

print_out=print
def get_derivatives_from_energy(the_di,t, print_info=False):
    # FD formula (polarisabilities are positive by convention)",
    plus, zero, minus= ( the_di[x] for x in  [t+'_pos', 'zero',t+'_neg'] )
    delta = (plus['delta'] - zero['delta'] )
    if print_info:
        print_out(f"\nE_plus={plus['E']} , E0={zero['E']}, E_minus={minus['E']} for tag={{t}} and delta={delta}")
    epl,eze,emi=[float(x['E']) for x in [plus,zero,minus]]
    return (
        get_first_derivative(epl,emi,delta) ,
        -get_second_derivative(epl,emi,eze,delta) 
    )
def get_derivatives_from_moments(the_di,t):
    """ Could also get second derivative """
    cur_pol=dict()
    try:
        plus, zero= ( the_di[x] for x in  [t+'_pos', 'zero'] )
    except: raise Exception(f" {the_di.keys()}")
    delta = (plus['delta'] - zero['delta'] )
    if any(x['mom'] is None for x in [plus, zero] ): return None
    for k in set(list(plus['mom'].keys())+list(zero['mom'].keys()) ):
        assert k in plus['mom'].keys() and k in zero['mom'].keys()
        cur_pol.update({ 
            k: np.array( -get_first_derivative(plus['mom'][k],zero['mom'][k], delta=delta/2) )
        })
    return cur_pol

def get_moments(dummy): pass
def make_entries_for_fields(pot, delta, do_mom):
    """ get_moments returns info di which may have additional keys for energy (starting with E_, e.g. Mp2)
    Transform info_di so that it return a dictionary with keys for method (e.g. 'HF','MP2', etc) each value containing a moments and and energy (currently only leading 'mom' available)
    """
    info_di, wfn, mom = get_moments( pot, delta, do_mom=do_mom)
    # get all methods (follow with _ after either E or mom)
    methods=[ k.split('_')[1] for k in info_di.keys() if '_' in k ]
    out_di=dict()
    for m in methods:
        out_di[m]=dict(
            E= info_di.get(f"E_{m}", None),
            mom= info_di.get(f"mom_{m}", None),
            delta= delta,
        )
    out_di['MAIN']=dict(
        E= info_di.get(f"E", None),
        mom= info_di.get(f"mom", None),
        delta= delta,
    )
    return out_di, wfn

def perform_df_calc(the_di):
    """Compute finite-difference derivatives and return dictionaries.

    Returns a dict with keys:
        - 'pol_eng': pol_through_energy
        - 'pol_mom': pol_through_moments (empty dict if gradients unavailable)
        - 'mom': mom_through_energy
    """

    pol_through_energy, mom_through_energy, pol_through_moments = dict(), dict(), dict() 

    zero_di=the_di['zero']
    assert len(zero_di)==1, f"Expected only one entry for zero field, got {len(zero_di)} with entries {zero_di}"
    zero_di=zero_di[0]
    for t in [ x for x in the_di.keys() if x.lower()!='zero']: # Tensor elements!
        # energy-based derivatives, 
        entries= the_di[t]
        # each dict of mom, E, delta 
        # unzip a sequence of (E, mom, delta) tuples produced per entry
        energies, moments, deltas = zip(*((e.get('E'), e.get('mom'), e.get('delta')) for e in entries))
        assert len(energies)==2, f"Expected 2 energy entries for plus and minus field, got {len(energies)} for tag {t} with entries {entries}"
        assert deltas[0]==-deltas[1], f"Expected deltas to be symmetric for plus and minus field, got {deltas} for tag {t} with entries {entries}"
        pol_through_energy[f"{t}@{t}"] = -get_second_derivative( energies[0], energies[1], zero_di.get('E'), deltas[0]-zero_di.get('delta') ) # negative due to convention
        mom_through_energy[f"{t}"] = get_first_derivative( energies[0], energies[1], deltas[0]-deltas[1] )

        found=([ (x,delta) for x,delta in zip(moments, deltas) if not x is None])
        # If one do two simple FF if two do both sides (removes square term)
        if len(found)==0:
            pol_through_moments[f"{t}@{t}"] = None
        elif len(found)==1:
            pol_through_moments[f"{t}@{t}"] = get_first_derivative(found[0][0], zero_di['mom'], found[0][1]-zero_di['delta'])
        elif len(found)==2:
            pol_through_moments[f"{t}@{t}"] = get_first_derivative(found[0][0], found[1][0], found[0][1]-found[1][1])

        # optionally get polarisabilities from moments if available

    return dict(
        pol_eng= pol_through_energy,
        pol_mom= pol_through_moments,
        mom= mom_through_energy
    )



@val_call
def generate_pol_ff_job_lines(job:psi4_job)->Tuple[lis,lis]:
    """Generate job lines to compute polarisability tensors and store them in a yaml file.

    Args:
        job (psi4_job): The psi4 job object.
    """
    lines=[]
    func_lines=[]
    #    ([] if not job.method.method.lower().startswith('cc') else [ ind+f"cc_convergence={opts.cc_convergence},"]) 
    lines+=[
        f"# For HF reference calculation as I understand. Should be probably a bit tighter for finite field (possibly its even faster do to better guess?)",
        f"e_convergence_val={job.job_settings.e_convergence}",
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
        f"delta_dip={job.job_settings.delta_dip:E}",
        f"delta_qad={job.job_settings.delta_qad:E}",
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
# 
    

    # make a list of dictionaries with ( <pot_tag>: ( delta=..., method1=(E:...,mom:....))'field': field_tag, 'deltas': [dots])
    lines+= [
        *(inspect.getsource(make_entries_for_fields).split('\n')),
        f"# Fields plus their finite difference values ",
        f"fields=['zero']+tensor_elements",
        f"deltas=[ [0] ]+[ [1*(delta_dip if len(t)==1 else delta_qad),-1*(delta_dip if len(t)==1 else delta_qad)] for t in tensor_elements]",
        f"save_fchk=[ [True] ]+[ (False,False) for _ in tensor_elements]",
        f"number_of_calculations= 1 + 2*len(tensor_elements) # zero plus positive and negative for each element",
        #
        f"# Evaluate wfn in external potential",
        f"{get_var('ff_reg')}=dict()",
        f"cnt=0",
        f"for i, (pot_tag,deltas, save) in enumerate(zip(fields,deltas,save_fchk)):",
        ind+f"for (delta,sv) in zip(deltas,save):",
        (ind+ind+f"do_mom=True if delta>0 else False  # only do moments for positive field (negative is for energy derviation) "
         if (job.job_settings.do_moments and job.method.can_do_gradient)
         else ind+ind+f"do_mom=False  #no moments request by main script"),
        ind+ind+f"info_di,wfn = {make_entries_for_fields.__name__}( {get_var('pot_registry')}[pot_tag] , delta, do_mom=do_mom )",
        ind+ind+f"if not pot_tag in {get_var('ff_reg')}.keys(): {get_var('ff_reg')}[pot_tag]=list()",
        ind+ind+f"{get_var('ff_reg')}[pot_tag].append(info_di)",
        ind+ind+f"cnt+=1",
        ind+ind+f"print_out(f\"CUSTOM MARKER BY USER: Finished moment '{{pot_tag}}' with delta={{delta}} ({{cnt}} out of {{number_of_calculations}}, {{(i+1)/number_of_calculations*100}}%)\")",
        ind+ind+f"if sv: # save fchk only for zero field and positive field (negative is for energy derviation) ",
        ind+ind+ind+f"try:",
        ind+ind+ind+ind+f"fchk_file='{job.jobname}.fchk'",
        ind+ind+ind+ind+f"psi4.fchk(wfn, fchk_file)",
        ind+ind+ind+ind+f"if not os.path.isfile(fchk_file): raise Exception(f\"fchk file was not created\")",
        ind+ind+ind+ind+f"print_out(f\"CUSTOM MARKER: Saved fchk file for '{{pot_tag}}' with delta={{delta}} as '{{fchk_file}}'\")",
        ind+ind+ind+f"except Exception as e:",
        ind+ind+ind+ind+f"print_out(f\"Could not save fchk file for '{{pot_tag}}' with delta={{delta}} due to {{e}}.\")",   
        ind+ind+ind+ind+f"fchk_file=None",
        * [ ind+ind+ind+x for x in store(job, ('files', f'final_fchk'), 'fchk_file')],# # save None if saving failed",
    ]

    # Now reshape the dictionary from entries of dictionaries to keys of 'MAIN','HF','MP2', etc each containing dictionary with the potentials tags and deltas and 'E' and 'mom' (if available)
    lines+= [
        f"# Reshape the dictionary to have method keys (e.g. 'HF','MP2', etc) each containing a dictionary with potential tags and their info (E and mom)",
        f"raw_di=dict()",
        #f"raise Exception({get_var('ff_reg')}) # DEBUG",
        f"for pot_tag, pot_infos in {get_var('ff_reg')}.items():",
        ind+f"for info_di in pot_infos:",
        ind+ind+f"for method, method_info in info_di.items():",
        ind+ind+ind+f"if method not in raw_di.keys(): raw_di[method] = dict()",
        # there might be different delta so pot_tag carries a list of dictionaries
        ind+ind+ind+f"if not pot_tag in raw_di[method].keys(): raw_di[method][pot_tag] = list()",
        ind+ind+ind+f"raw_di[method][pot_tag]+=[ dict( (k,v) for k,v in method_info.items() if k in ['E','mom','delta'] ) ]",
    ]

    # Compute the finite difference

    # Now perform the DF calculation, each method tag will results in a entry in pol_reg with the polarisabilities through energy and through moments (if available)
    var_mem.update({f"pol_reg":"pol_reg"})
    lines+=[
        f"",
        f"# Compute polarisabilities through finite differences",
        f"# Each method (e.g. 'HF','MP2', etc) will have an entry in {get_var('pol_reg')} with the polarisabilities through energy and through moments (if available)",
        f"",
        f"#FORMULAS",
        *inspect.getsource(get_second_derivative).split('\n'),
        *inspect.getsource(get_first_derivative).split('\n'),
        *inspect.getsource(get_derivatives_from_energy).split('\n'),
        *( [] if not job.method.can_do_gradient else 
            inspect.getsource(get_derivatives_from_moments).split('\n')
        ),
        *inspect.getsource(perform_df_calc).split('\n'),
    ]+[
        f"#BODY",
        f"{get_var('pol_reg')}=dict()",
        f"for method in raw_di.keys():",
        # treat the zero field case as well, to get the energy and moments at zero field (which may be different from the reference if e.g. MP2 correlation energy is only available as correction to HF)
        # ind+f"zero = raw_di[method].pop('zero') # remove the zero field entry from raw_di and save it separately (we need it for the derivatives and it should not be in the final pol_reg)",
        ind+f"{get_var('pol_reg')}[method] = {perform_df_calc.__name__}(raw_di[method])",
    ]



    path=('results','properties')
    lines+=[f"#Store object",  
            *store(job,path+('expansion_center',), f"[ float(x) for x in {get_var('center_of_charge')}]"),
    ]
    # Now store the different energy entries. MAIN_ is just mapped on nothing the remained is used a prefix (HF_, MP2_, etc) for the method. The pol_reg contains the polarisabilities through energy and through moments (if available)
    lines+=[
        f"for method, values in {get_var('pol_reg')}.items():",
        ind+f"method_prefix = method+'_' if method != 'MAIN' else '' # MAIN is the reference method, so we don't add a prefix to it",
        ind+f"for key, val in values.items():",
        # map keys onto MolPol_FinFie_through_eng and MolPol_FinFie_through_dens, and for the moments MolMom_FinFie_through_eng and MolMom_FinFie_through_dens (if available)
        ind+ind+f"if key.startswith('pol_eng'): ",
        ind+ind+ind+f"store_key = f\"{{method_prefix}}MolPol_FinFie_through_eng\"",
        ind+ind+f"elif key.startswith('pol_mom'): ",
        ind+ind+ind+f"store_key = f\"{{method_prefix}}MolPol_FinFie_through_dens\"",
        ind+ind+f"elif key.startswith('mom'): ",
        ind+ind+ind+f"store_key = f\"{{method_prefix}}MolMom_FinFie_through_eng\"",
        *[ ind+ind+x for x in store(job,(*path,"__store_key__"), '( val if any( x is not None for x in val.values() ) else None )')],
    ]

    if do_single_point_mom:
        lines+=[
            f"# Optionally store the moments as well (only for single point with gradient available, as we only get the moments for the positive field, which are not suitable for finite difference derivatives but may be interesting on their own)",
            *(store(job,path+('MolMom_from_dens',), f"raw_di['MAIN']['zero'][0]['mom']")),
        ]



    return lines, func_lines



# set_options( {'perturb_dipole': make_c2v_vec(delta),'reference':'rhf' })
# set {
# OPDM_RELAX True
# }
