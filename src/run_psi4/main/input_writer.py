from . import *

def set_hardware(memory='4GB', num_threads=4, storage_var=None):
    hardware_settings=dict(
        memory=p4ut.embed_value(memory),
        num_threads=p4ut.embed_value(num_threads),
    )
    content=[]
    content+=[f"# Hardware settings "]
    hardware_lines=[ 4*" "+f"{k}={v}," for k, v in hardware_settings.items() ]
    content+=["hardware_settings=dict("]+ hardware_lines +[ ")"]
    content+=[f"set_{k}(hardware_settings['{k}'])" for k ,v in hardware_settings.items()]
    
    if not isinstance(storage_var, type(None)):
        content+=[ f"{storage_var}['about_the_run'].update(hardware_settings)"]
    content='\n'.join([ f"{line}" for line in content])
    return content

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
def molecule_string(atom_types, coordinates, units=None, total_charge=0, multiplicity=1):

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

    return '\n'.join(content)

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
def get_polarizability(method, geom=None, storage_var=None):
    """ Given a geometry calculate the Polarisabilities available in Psi4
    [ ] Provide starting guess (idk if that will be that much faster though)
    """

    # Get method information
    func=method['functional']
    basis_set=method['basis_set']
    # Get geometry information 
    molecule_name=geom['molecule_name']
    atom_types=geom['atom_types']
    coordinates=geom['coordinates']

    def calculate_center_of_mass(atom_types, coordinates):
        """ 
        Calculate center of mass as
        x_com = sum ( w_i * x_i )    with x_i beeing coordinates and w_i weights: w_i = Z_i/( sum Z_j )
        """
        atom_numbers=[]
        for at_ty in atom_types:
            if isinstance(at_ty, int):
                atom_numbers.append(at_ty)
            elif isinstance(at_ty, str):
                assert len(at_ty)<3 and len(at_ty)>0
                Z=element_symbol_to_nuclear_charge(at_ty)
                atom_numbers.append(Z)
        atom_numbers=np.array(atom_numbers)
        com=np.sum(atom_numbers[:,None]*coordinates, axis=0)/np.sum(atom_numbers)
        return com


    COM, units_of_origin=[ [ float(x) for x in calculate_center_of_mass(atom_types, coordinates)], 'ANGSTROM' ] # recast that for easy printout
        

    property_calculation=[
    "# perform polarizability calculation",
    "props = ['DIPOLE_POLARIZABILITIES','QUADRUPOLE_POLARIZABILITIES']",
    #"set_options({",
    #"   \'basis\' : basis_set,",
    #"   \'reference\' : 'rks',",
    #"})",
    "print(method, monomer, props)",
    "E, wfn = properties(method, molecule=monomer, properties= props, return_wfn=True)",
    ]

    # Get the polarisabilities
    def recover_polarisabilities(origin, type_of_origin=None, units_of_origin=None):
        recover_polarisabilities=f"""
        #### Recover the polarisabilities saved in the global variables
        # Manipulate the 'variables' variable where polarisabilities are stored
        the_variables=variables()
        key_per_rank={{
            'dipole':'DIPOLE POLARIZABILITY',
            'quadrupole':'QUADRUPOLE POLARIZABILITY'
        }}
        # Construct all expected subindices
        def make_product(the_list):
        \tnew_list=[]
        \tfor x in the_list:
            \t\tfor y in the_list:
            \t\t\tnew_list.append(x+y)
        \treturn new_list
        expected_subindices={{
            'dipole':make_product(['x','y','z']),
            'quadrupole':make_product(['xx','yy','zz','xz','xy','yz']),
        }}

        polarisability_tensors={{}}
        zero_thres=1.e-12
        for rank, rank_key in key_per_rank.items():
        \ttensor={{}}
        \tavailable_keys=[x for x in the_variables.keys() if x.startswith(rank_key)]
        \tcovered_keys=dict([(x,False) for x in available_keys ])
        \tfor index in expected_subindices[rank]:
            \t\tfull_index=' '.join([rank_key,index.upper()])
            \t\tassert full_index in available_keys, [index, available_keys]
            \t\t# Mark the index as covered
            \t\tassert not covered_keys[full_index], f"Double read in of {{full_index}}?"
            \t\tcovered_keys[full_index]=True
            \t\t# Safe the index and value into the tensor
            \t\tvalue=the_variables[full_index]
            \t\tvalue=0 if abs(value)<zero_thres else value
            \t\ttensor.update({{index:value}})
        \tpolarisability_tensors.update({{rank:tensor}})
        \tassert all(covered_keys.values())

        # Drop to file
        polarisability_tensors.update({{ 'origin': {{'coordinate':{COM}, 'type_of_origin':\'{type_of_origin}\', 'units_of_origin':\'{units_of_origin}\' }}}})
        polarizability_file='polarisability_tensors.yaml'
        with open(polarizability_file,'w') as wr:
          \tyaml.dump(polarisability_tensors, wr)
        {storage_var}['files'].update({{'polarizability_tensors':polarizability_file}})
        """
        ind=4*' '
        recover_polarisabilities=[x.strip(' ') for x in recover_polarisabilities.split('\n')]
        return recover_polarisabilities
        # Sorting according to stones convention (See Cartesian–Spherical Conversion Tables@Stone's Theory of Intermolecular Interactions book)

    define_polarisability_origin=[x.strip() for x in f"""
    # Define center of mass (in angstrom!) and assert that it has been read in correctly!
    core.set_global_option("PROPERTIES_ORIGIN", {COM})
    assert {COM}==get_option('None','PROPERTIES_ORIGIN'), 'Check Center of mass'
    """.split('\n')]
    # Pierce file together
    lines=(        []
        + define_polarisability_origin
        + property_calculation
        + recover_polarisabilities(COM, type_of_origin='COM',units_of_origin=units_of_origin)
    )
    return lines

def get_gs_energy(func, basis, molecule_name='monomer', uks=False):
    E_var='E_neut'
    wfn_var='wfn_neut'
    if uks:
        ref='uks'
    else:
        ref='rks'
    base_string= [ x.strip() for x in f"""
    set reference {ref}
    set_options({{
        'basis': {basis},
        'reference': 'rks',
    }})
    {E_var},{wfn_var} = energy({func}, molecule={molecule_name}, return_wfn=True)
    """.split('\n') ]
    return base_string, E_var, wfn_var

def get_ionization_energy(func_var, basis_set_var, molecule_name='monomer', uks=True):
    E_var='E_ion'
    wfn_var='wfn_ion'
    if uks:
        ref='uks'
    else:
        ref='rks'
    molecule_ion_name=molecule_name+"_ion"
    base_string= [ x.strip() for x in f"""
    {molecule_ion_name}=Molecule.from_string({molecule_name}_string)
    {molecule_ion_name}.set_molecular_charge({molecule_name}.molecular_charge()+1)
    {molecule_ion_name}.set_multiplicity({molecule_name}.multiplicity()+1)
    set_options({{
        'basis': {basis_set_var},
        'reference': 'uks',
    }})
    {E_var},{wfn_var} = energy({func_var}, molecule={molecule_ion_name}, return_wfn=True)
    """.split('\n') ]
    return base_string, E_var, wfn_var

def get_grac_shift(func_var, basis_set_var, molecule_name='monomer', charge=0, spin=1):
    neut_string, E_neut, wfn_neut = get_gs_energy(func_var, basis_set_var, molecule_name=molecule_name)
    ion_string, E_ion, wfn_ion    = get_ionization_energy(func_var, basis_set_var, molecule_name=molecule_name)

    # Compute the difference -> grac_shift
    ac_var=f"AC_shift"
    string=[ x.strip() for x in f"""
    E_homo = float({wfn_neut}.epsilon_a_subset('AO', 'ALL').np[{wfn_neut}.nalpha()-1])
    IP = float({E_ion}) - float({E_neut})
    {ac_var} = IP + E_homo
    """.split('\n') ]

    store_string=[
        f"storage['results'].update({{ 'HOMO': E_homo, 'IP': IP , 'AC_SHIFT': {ac_var} }})",
    ]

    grac_string=neut_string+ion_string+string+store_string
    return grac_string, ac_var
def set_grac_shift(ac_var:str|float):
    return [x.strip() for x in f"""
        set_options({{'dft_grac_shift': {ac_var}}})
        the_grac_shift=get_option('SCF','dft_grac_shift')
        print_out(f\"GRAC SHIFT SET TO: {{the_grac_shift}}\")
        assert the_grac_shift>0.0, f\"The grac shift is not larger zero!\"
        """.split('\n')
    ]
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