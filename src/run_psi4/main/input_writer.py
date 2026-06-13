from . import *


def molecule_input(atom_types, coordinates, units=None, total_charge=0, multiplicity=1, name='monomer'):

    if isinstance(units, type(None)):
        units={'LENGTH':'BOHR'}
    content=[]

    # Comment (makes it easier to copy this block into xyz
    content+=[f"# Number of nuclei: {len(coordinates)}"]

    # Charge and multiplicity
    content+=[f"{total_charge} {multiplicity} # Charge and Multiplicity"]
    
    # Make geometry
    geometry=[]
    for at_ty, coord in zip(atom_types, coordinates):
        assert len(at_ty)<3 and at_ty!=''
        coord=[ f"{c:.8f}" for c in coord]
        coord=[ f"{c:>15}" for c in coord]
        geometry+=[ f"{at_ty.capitalize():<2} {' '.join(coord)}"]
    content+=geometry


    content+=[
        'Symmetry C1',
        'nocom',
        'noreorient',
        f"units {''.join(list(units.values()))}"
    ]

    box=f"molecule {name} {{\n__REP__}}\n"
    spacer=4* ' '
    content=''.join([ f"{spacer}{line}\n" for line in content])
    string=box.replace('__REP__', content)
    return string

def calc_string(dft_functional, basis, molecule, reference='rks', en_var='E', wfn_var='wfn', marker='NO_MARKER', do_timing=True,
                #save_wfn=False, jobname=None, tag=None
    ):
    
    #if save_wfn:
    #    assert not isinstance(jobname,type(None)), f"Save wave function requested but jobname not provided!"

    freeze_core=True
    content=[]
    if do_timing:
        content+=[
            f"start=time()",
        ]
    content+=[
        f"set reference {reference}",
        f"set basis {basis}",
        f"set freeze_core {freeze_core}",
    ]
    content+=[ f"{en_var}, {wfn_var} = energy(\"{dft_functional}\", molecule={molecule}, return_wfn=True)" ]

    if do_timing:
        identifier=f"TIMING {marker:<15}"
        content+=[
            f'# TIMING,',
            f'duration=time()-start' ,
            f"psi4.print_out(  \'\\n\'.join([",
                f"f\"{identifier} {'[in seconds]':<15} : {{duration}}\",",
                f"f\"{identifier} {'[in H:M:S]':<15} : {{int(duration//3600)}}-{{int(duration%3600//60):02d}}-{{duration%60}}\",",
            '])  )'
        ]
    #if save_wfn:
    #    content+=save_wfn_string(wfn_var=wfn_var, jobname=jobname, tag=tag)


    string=''.join([f"{x}\n"  for x in content])
    return string

def grac_shift(E_ion, E_neut, wfn_neut, verbose=True):
    ac_var=f"AC_shift"
    homo_var=f"E_homo"
    
    string=[     f"{homo_var} = {wfn_neut}.epsilon_a_subset('AO', 'ALL').np[{wfn_neut}.nalpha()-1]"
        ,f"IP = {E_ion} - {E_neut}"
        ,f"{ac_var} = IP + {homo_var}"
        ,f"set dft_grac_shift {ac_var}"
    ]

    if verbose:
        ind=4*' '
        information_string=[f"Information about GRAC calculation",
            f"IP {{IP}}",
            f"E_homo {{{homo_var}}}",
            f"GRAC_shift {{{ac_var}}}",
        ]
        information_string='\n'.join(
            [f"\npsi4.print_out(  \'\\n\'.join([",]+
            [ind+'\"'+x+'\",' for x in information_string]+
            ["])  )",]
        )
        string+=[ information_string ]
    string='\n'.join(string)
    return string

def psi4_settings(output_file):
    local_content=['# Non-hardware related psi4 settings']
    local_content.append(f"set_output_file(\"{output_file}\")")
    return '\n'.join(local_content)

def save_wfn_string(wfn_var, jobname, tag, save_fchk=True, save_npy=False):
    """ Return strings that saves wave function corresponding to provided variable name """
    local_content=[f"# Save wave function to file"]
    ### Save wave function
    if isinstance(tag, type(None)):
       tag=''
    else:
        tag=f"_{tag}" 

    if save_npy:
        wfn_file=f"{jobname}{tag}.wfn.npy"
        local_content+=[f"{wfn_var}.to_file(\"{wfn_file}\")"]
    else:
        wfn_file=None

    if save_fchk:
        fchk_file=f"{jobname}{tag}.fchk"
        local_content+=[f"fchk({wfn_var},\"{fchk_file}\")"]
    else:
        fchk_file=None

    return local_content, wfn_file, fchk_file
    return lines

@val_call
def get_gs_energy(func, basis, molecule_name='monomer', uks=False) ->Tuple[List[str],str,str]:
    E_var='E_neut'
    wfn_var='wfn_neut'
    if uks:
        ref='uks'
    else:
        ref='rks'
    base_string= [
    f"set_options({{",
    # ind+f"\'basis\': {basis},",
    ind+f"'reference': '{ref}',",
    f"}})",
    ]
    base_string+=[
    f"{E_var},{wfn_var} = energy({func}, molecule={molecule_name}, return_wfn=True)"
    ]
    return base_string, E_var, wfn_var

@val_call
def get_ionization_energy(func_var, basis_set_var, molecule_name='monomer', molecule_string='molecule_string', uks=True,
                          conv_crit=None) -> Tuple[List[str],str,str]:
    """ 
    conv_crit: currently only None provide """
    E_var='E_ion'
    wfn_var='wfn_ion'
    if uks:
        ref='uks'
    else:
        ref='rks'

    molecule_ion_name=molecule_name+"_ion"
    base_string= [ x.strip() for x in f"""
    {molecule_ion_name}=Molecule.from_string({molecule_string})
    assert {molecule_name}.multiplicity()==1, \"Expected multiplicity 1 for ion, got {molecule_name}.multiplicity (need to handle other ones)\"
    {molecule_ion_name}.set_molecular_charge({molecule_name}.molecular_charge()+1)
    {molecule_ion_name}.set_multiplicity({molecule_name}.multiplicity()+1)
    set_options({{
        'reference': 'uks',
    }})
    {E_var},{wfn_var} = energy({func_var}, molecule={molecule_ion_name}, return_wfn=True)
    """.split('\n')
    ]
    if conv_crit is None:
        pass
    elif isinstance(conv_crit, str) and conv_crit=='loose':
        get_conv=[
            "old_conv=dict( (k, get_global_option(k))  for k in [ f\"{x}_convergence\" for x in ['e','d','r'] ])",
            "ion_conv=dict( (k, 1.e-6) for k in old_conv.keys() ) # usually ion is calculated for grac and does not need to be tight",
            'set_options(ion_conv)'
        ]
        reset_conv=[
            f"set_options(old_conv)"
        ]
        base_string=[ *get_conv,*base_string,*reset_conv ]
    else:
        raise NotImplementedError(f"Noth convergence control implemented for {get_ionization_energy} and value {conv_crit}")
    
    return base_string, E_var, wfn_var

    # if verbose:
    #     ind=4*' '
    #     information_string=[f"Information about GRAC calculation",
    #         f"IP {{IP}}",
    #         f"E_homo {{{homo_var}}}",
    #         f"GRAC_shift {{{ac_var}}}",
    #     ]
    #     information_string='\n'.join(
    #         [f"\npsi4.print_out(  \'\\n\'.join([",]+
    #         [ind+'\"'+x+'\",' for x in information_string]+
    #         ["])  )",]
    #     )
    #     string+=[ information_string ]
    # string='\n'.join(string)
    # return string

def get_ip_tuned_omega(func, basis, molecule_name='monomer', charge=0, spin=1, storage_var=None):
    """IP tuning requires calcultion of ionization energy [do this once in this case]
    then we alter omega so that we obtain a homo value that equals the negative ionization energy
    """
    assert not isinstance(storage_var, type(None)), f"Provide storage var!"
    neut_string, E_neut, wfn_neut = get_gs_energy(func, basis, molecule_name=molecule_name)
    ion_string, E_ion, wfn_ion    = get_ionization_energy(func, basis, molecule_name=molecule_name)

    # Compute the difference -> grac_shift
    ac_var=f"AC_shift"
    get_IP_string=[ x.strip() for x in f"""
    IP = {E_ion} - {E_neut}
    print_out(f"The calculate IP is: {{IP}}")
    {storage_var}['results'].update({{'IP':IP}})
    """.split('\n') ]


    optimize_omega_func=f"""
    def get_ip_diff(molecule, func, basis, omega, IP):
        set_options({{
            'reference' : 'rks',
            'DFT_omega' : omega,
            'basis'     : basis,
        }})
        E_omega, wfn_omega = energy(func, molecule=molecule, return_wfn=True)
        E_homo = wfn_omega.epsilon_a_subset('AO', 'ALL').np[wfn_omega.nalpha()-1]
        diff=E_homo+IP
        return diff
    """
    optimization_algorithm=f"""
    def bisection(x_vals, y_vals):
        x_lw, x_md, x_up =x_vals
        new_x_lw, new_x_md, new_x_up =x_vals
        y_lw, y_md, y_up = y_vals

            
        
        if y_lw*y_md<0:
            new_x_up=x_md
            new_x_lw=x_lw
        else:
            new_x_up=x_up
            new_x_lw=x_md
        new_x_md=(new_x_up+new_x_lw)/2.
        return  [new_x_lw, new_x_md, new_x_up]
    def return_diff_values(omega_range, old_omega_range, old_omega_vals):
        if isinstance(old_omega_range, type(None)):
            new_vals = [ get_ip_diff({molecule_name}, method, basis_set, om, IP) for om in omega_range]
        else:
            mapper=dict( [ (k,v) for k,v in zip(old_omega_range, old_omega_vals)])
            new_vals=[]
            assert sum([ x in mapper.keys() for x in omega_range ]) == 2, f"{{mapper}} {{omega_range}}"
            for om in omega_range:
                if om in mapper.keys():
                    new_vals.append(mapper[om])
                else:
                    new_vals.append(get_ip_diff({molecule_name}, method, basis_set, om, IP))
        return new_vals

    omegas_ini=[ 0.1 ,  .5 , .9 ]
    omegas=omegas_ini
    old_omegas=None
    old_diff_vals=None
    chosen_omega=None
    for i in range(10):
        thres=0.005
        if np.abs(omegas[0]-omegas[2])<thres:
            chosen_omega=omegas[1]
            break
        else:
            diff_vals=return_diff_values(omegas, old_omegas, old_diff_vals)
            if diff_vals[0]*diff_vals[1]>=0. :
                test_vals=[0.05*i for i in range(1,20)]
                # raise Exception( test_vals, [ float(get_ip_diff({molecule_name}, functional, basis_set, om, IP)) for om in test_vals ])
            old_omegas=omegas
            old_diff_vals=diff_vals
            print(omegas, diff_vals)
            omegas=bisection(omegas, diff_vals)
            print('new_omegas',omegas)
    assert not isinstance(chosen_omega, type(None))
    {storage_var}['results'].update({{"chosen_omega":chosen_omega}})
    """

    def splitter(string):
        lines=string.split('\n')
        first_line=None
        for i in range(len(lines)):
            if lines[i].strip()=='':
                continue
            else:
                first_line=i
                break
        assert not isinstance(first_line, type(None))
        indent=( len(lines[first_line])-len(lines[first_line].lstrip(' ')) )*' '
        return [l.replace(indent,'',1) for l in lines[first_line:]]
    optimize_omega_func=splitter(optimize_omega_func)
    optimization_algorithm=[f"start_omega_optimization=time.time()",]+splitter(optimization_algorithm) \
        + [f"{storage_var}['about_the_run']['timings'].update({{'omega_optimization':time.time()-start_omega_optimization}})"]

    # {ac_var} = IP + E_homo

    get_omega_str=[f"start_get_IP=time.time()"]+neut_string+ion_string+get_IP_string+[f"{storage_var}['about_the_run']['timings'].update({{'get_IP':time.time()-start_get_IP}})"]+optimize_omega_func+optimization_algorithm
    omega_var=f"chosen_omega"

    return get_omega_str, omega_var

def set_omega(omega_var):
    return [f"set DFT_omega {omega_var}",
        f"print_out(f\"IP-TUNED OMEGA SET TO: {{{omega_var}}}\")"]