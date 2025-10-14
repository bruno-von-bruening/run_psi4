
import os, sys, re, time
import numpy as np


from .input_writer import *

def psi4_exit_ana(output_file):
    with open(output_file, 'r') as rd:
        lines=rd.readlines()

    normal_exit_line_fl=bool(re.match(r'\*\*\* Psi4 exiting successfully[.]*^', lines[-1]))

    if normal_exit_line_fl:
        exit_status='normal'
    else:
        exit_status='failed'
    return {
        'exit_status': exit_status
    }
                    
def complete_calc(
    atom_types, coordinates, total_charge=0, multiplicity=1, # structure
    dft_functional='PBE0', basis_set='aug-cc-pVTZ', do_GRAC=True, # calculation
    jobname='test', units=None,
    hardware_settings={'memory':'4GB', 'num_threads':4}
):
    params=dict(
        save_standard_wfn=True,
    )
    save_npy=False
    content=[]

    # Comment
    the_time=f"{time.strftime('%Y-%M-%D_%H:%M:%S', time.localtime())} (Zone: {time.tzname})"
    comment_string=f"# Psi4 calculation for job {jobname} input file\n# generated at {the_time}"
    content+=[comment_string, '']

    content+=[f"from time import time",'']

    #Set the hardware
    hardware_string=set_hardware(**hardware_settings)
    content+=[hardware_string,'']

    # Set psi4 settings
    output_file=f"{jobname}.psi4out"
    psi4_settings_string=psi4_settings(output_file=output_file)
    content+=[psi4_settings_string,'']

    ## Set the molecule data
    dat_monomer={
        'name':'monomer',
        'total_charge':total_charge,
        'multiplicity':multiplicity,
    }
    string=molecule_input(atom_types, coordinates, **dat_monomer, units=units)
    content+=[string]
    

    if do_GRAC:
        # Shift charge and multiplicity
        charge_ion=total_charge+1
        if multiplicity!=1:
            raise Exception(f"Can only handle multiplicty of 1 for GRAC at the moment!")
        else:
            multiplicity_ion=2

        dat_monomer_ion={
            'name':         'monomer_ion',
            'total_charge':       charge_ion,
            'multiplicity': multiplicity_ion,
        }
        string_ion=molecule_input(atom_types, coordinates, **dat_monomer_ion, units=units)
        content+=[string_ion]

    #### Calculations
    reference='rks' if multiplicity==1 else 'uks'
    calc_neut={
        'molecule':'monomer',
        'reference':reference,
        'en_var':'E_neut',
        'wfn_var':'wfn_neut'
    }
    string_calc_neut = calc_string(dft_functional, basis_set, **calc_neut, marker='SCF_NEUTRAL', reference=reference)
    save_calc_neut, wfn_neut_file, fchk_neut_file= save_wfn_string(calc_neut['wfn_var'], jobname=jobname, tag='neut', save_npy=save_npy) 
    content+=[f"#Calculate neutral molecule", string_calc_neut, *save_calc_neut]


    if do_GRAC:
        calc_ion={
            'molecule':'monomer_ion',
            'reference':'uks',
            'en_var':'E_ion',
            'wfn_var':'wfn_ion'
        }
        string_calc_ion = calc_string(dft_functional, basis_set, **calc_ion, marker='SCF_ION')
        content+=[f"#Calculate ion", string_calc_ion]

        # Set acs
        ac_string = grac_shift( E_ion=calc_ion['en_var'], E_neut=calc_neut['en_var'], wfn_neut=calc_neut['wfn_var'])

        calc_grac={
            'molecule':'monomer',
            'reference':'rks',
            'en_var':'E_grac',
            'wfn_var':'wfn_grac'
        }
        string_calc_grac = calc_string(dft_functional, basis_set, **calc_grac, marker='SCF_GRAC')
        save_calc_grac, wfn_grac_file, fchk_grac_file= save_wfn_string(wfn_var=calc_grac['wfn_var'], jobname=jobname, tag='grac', save_npy=save_npy) 
        content+=[f"# GRAC SHIFT CALCULATION\n"+ac_string+f"\n# GRAC CORRECTED SCF\n"+string_calc_grac, *save_calc_grac]
    


    # Remove temporary files
    content+=['clean()']

    ### Print file
    psi4inp_path=f"{jobname}.psi4inp"
    content='\n'.join(content)
    with open(psi4inp_path,'w') as wr:
        wr.write(content)

    psi4_files= {
        'psi4inp_file': os.path.realpath(psi4inp_path),
        'psi4out_file': os.path.realpath(output_file),
        #'wfn_neut_file': os.path.realpath(wfn_neut_file),
        'fchk_neut_file': os.path.realpath(fchk_neut_file),
    }
    if do_GRAC:
        psi4_files.update({
            #'wfn_grac_file': os.path.realpath(wfn_grac_file),
            'fchk_grac_file': os.path.realpath(fchk_grac_file),
        })
    psi4_files.update({
        'do_GRAC':do_GRAC,
    })
    return psi4_files

def get_single_point(method, geom=None, storage_var=None, jobname=None):
    
    method=method['functional']
    meth_fl=p4ut.get_method_flags(method)

    molecule_name=geom['molecule_name']
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
    lines=com+[
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

    return lines

def xyz_file_read(file_name):
    with open(file_name, 'r') as rd:
        lines=rd.readlines()
    assert len(lines)>2
    coordinates=[ x.split()[1:4] for x in lines[2:] if x.strip()!='' ]
    atom_types=[ x.split()[0] for x in lines[2:] if x.strip()!='' ]
    num_atoms=int(lines[0])
    assert len(coordinates)==num_atoms, f"{len(coordinates)} {num_atoms}"
    coordinates=np.array([ list(map(float, x)) for x in coordinates])

    return atom_types, coordinates

def xyz_to_string(file_name):
    """
    """
    atom_types, coordinates = xyz_file_read(file_name)
    molecule_str=molecule_string(atom_types=atom_types, coordinates=coordinates, units={'LENGTH':'ANGSTROM'})

    return molecule_str




def process_method(func, basis_set):
    # Check for grac
    func_comp=func.upper().split('-')
    grac_comp=[x for x in  func_comp if x.upper() in ['GRAC'] ]
    do_grac=False
    if len(grac_comp)>0:
        func='-'.join([x for x in func_comp if not x in grac_comp ])
        do_grac=True
    
    # Check for IP-tuned
    ip_tags=['IPTUNED','IPTUNING']
    func_comp=func.upper().split('-')
    ip_comp=[x for x in func_comp if x.upper() in ip_tags]
    do_ip_tune=False
    if len(ip_comp)>0:
        func='-'.join([x for x in func_comp if not x in ip_comp])
        do_ip_tune=True


    def splitter(string):
        split=[ x.strip() for x in string.split('\n') ]
        return '\n'.join(split)

    return func, basis_set, do_ip_tune, do_grac

def psi4_input_single(job,func, basis_set, geom, open_spin:bool, input_file_name, storage_file_name=None, jobname=None,
    num_threads=4, memory='4GB',
):
    """ Create an single psi4 input file """

    # Storage variable
    def make_storage(storage_var, storage_file):
        # Storate creation
        storage_file_var="storage_file"
        storage_lines=[
            f"# Store variables here (eg IP)",
            f"{storage_var}={{'files':{{}},  'results':{{}}, 'settings':{{}}, 'about_the_run':{{ 'timings': {{}} }} }}",
            f"{storage_file_var}=\'{storage_file}\'",
        ]
        return storage_lines, storage_file_var
    def drop_storage(storage_var, storage_file_var):
        # Storage dump
        storage_drop_lines=[
            f"# Time the complete run",
            f"{storage_var}['about_the_run']['timings'].update({{ 'total_run_time': time.time()-run_started_at,"+
                f" 'unaccounted_time': time.time()-run_started_at - sum(list({storage_var}['about_the_run']['timings'].values()))}})",
            f"# Drop the storage as yaml file",
            f"def drop_yaml(storage, storage_file):",
            f"    with open(storage_file, 'w') as wr:",
            f"        yaml.dump({storage_var},wr)",
            f"    print(f\"Dropped stored variables under {{storage_file}}\")",
            f"drop_yaml({storage_var}, {storage_file_var})",
        ]
        return storage_drop_lines
    
    imports=[
            f"# Imports",
            f"import time",
            f"import yaml",
            f"run_started_at=time.time()",
    ]

    # Care about the storage
    storage_var=f"storage"
    if isinstance(storage_file_name, type(None)):
        storage_file_name=f"psi4_calc_variable_storage.yaml"
    storage_lines, storage_file_var=make_storage(storage_var, storage_file_name)
    storage_drop_lines=drop_storage(storage_var, storage_file_var)

    # Creation of the molecul, jobname=Nonee
    molecule_lines=geom['molecule_lines']

    func, basis_set, do_ip_tuning, do_grac = process_method(func, basis_set)
    method=dict(
        functional=func,
        basis_set=basis_set,
        do_ip_tuning=do_ip_tuning,
        do_grac=do_grac,
    )

    func_var='method'
    basis_set_var='basis_set'
    definitions=[ x.strip() for x in f"""\
    # Description of the method (indicate correction in comments)
    {func_var}='{func}'
    {basis_set_var}='{basis_set}'
    # IP_TUNING: {do_ip_tuning}
    # GRAC correction: {do_grac}
    {storage_var}['settings'].update({{\'method\': method, 'basis_set': basis_set, \'grac_correction\':{do_grac}, \'ip_tuning\':{do_ip_tuning}}})"""
    .split('\n') ]


    # Set up the method
    if do_grac:
        if 'ac_shift' in geom.keys():
            if geom['ac_shift'] is None:
                grac_string, grac_var=get_grac_shift(func_var,basis_set_var)
            else:
                grac_var=geom['ac_shift']
                grac_string=[f"# Using user provided AC shift of {grac_var}"]
        else:
            grac_string, grac_var=get_grac_shift(func_var,basis_set_var)
        set_grac_str=set_grac_shift(grac_var)
    else:
        grac_string=[]
        set_grac_str=[]
    
    if do_ip_tuning:
        ip_tune_string, omega_var=get_ip_tuned_omega(func_var, basis_set_var, storage_var=storage_var)
        set_omega_str=set_omega(omega_var)
    else:
        ip_tune_string=[]
        set_omega_str=[]

    def make_method_definition_str(qm_protocol):

        method=qm_protocol['functional']
        meth_flags=p4ut.get_method_flags(method)

        do_ip_tuning=qm_protocol['do_ip_tuning']
        do_grac=qm_protocol['do_grac']

        the_options=dict(
            basis=basis_set_var,
        )
        # MP2 related keywords
        if meth_flags['mp2']:
            the_options.update( p4ut.make_mp2_options() )
        elif meth_flags['cc']:
            the_options.update( p4ut.make_cc_options() )
        elif meth_flags['hf']:
            the_options.update( p4ut.make_hf_options(open_spin=open_spin) )
        else:
            the_options.update( p4ut.make_dft_options(open_spin=open_spin) )
        
        if do_ip_tuning:
            the_options.update(dict(
                DFT_omega=f"{{{omega_var}}}"
            ))
        if do_grac:
            the_options.update(dict(
                DFT_GRAC_SHIFT=f"{grac_var}",
            ))
        
        # Print the options into the file
        indent=4*' '
        options_set=[ 
            "{indent}{k:<15} : {v:<15} ,".format(indent=indent, k='\''+k+'\'', v=v.strip('{}')) 
            for k,v in the_options.items() 
        ]
        # also print the options into the psi4 output 
        options_printout=[
            f"print_out(\"{indent}{k}={{}}\\n\".format({v}))" for k,v in the_options.items()
        ]

        # Create the whole section
        method_setup_lines=['# Setting the options for the calculation']
        method_setup_lines+=[ 'clean_options() # Set options to default (avoid bad surprises)', "set_options({"] + options_set + ["})"]
        info_lines=[f"# Information about the settings printed to output file", f"print_out(f\"INFORMATION about SET METHOD\\n\")"]
        info_lines+=options_printout
        
        return method_setup_lines+info_lines
    
    method_lines=make_method_definition_str(method)



    job_indication_lines=[
        10*"#",
        3*'#'+f" {job.upper()}",
        10*"#",
    ]
    hardware_lines=set_hardware(memory=memory,num_threads=num_threads, storage_var=storage_var).split('\n')
    current_time=time.strftime(f"%Y/%m/%d_%H:%m:%s_{time.tzname}", time.localtime())
    header=[
        f"# This is an auto-generated input file",
        f"# Generated by {__file__}",
        f"# time is {current_time}",
    ]
    content=( header +[''] + imports + [""] + storage_lines+[""] +hardware_lines+[""]+molecule_lines+[""]
        + definitions + ['']
        + ip_tune_string + grac_string
        +job_indication_lines+['']+method_lines+['']
    )
    # Get the job
    if job in ['polarizabilities']:
        get_polarizability_lines=get_polarizability(method, geom=geom, storage_var=storage_var)
        content+=['started_polarizabilities=time.time()'] + get_polarizability_lines + \
        [f"{storage_var}['about_the_run']['timings'].update({{'polarizabilities':time.time()-started_polarizabilities}})"]
    elif job in ['single_point','sp']:
        single_point_lines=get_single_point(method, geom=geom, storage_var=storage_var, jobname=jobname)
        content+=['# Starting a single point calculation', 'started_single_point=time.time()'] + single_point_lines 
        content+=[f"{storage_var}['about_the_run']['timings'].update({{'single_point':time.time()-started_single_point}})"]
    else:
        raise Exception(f"No job defined for {job}")

    content+=[""]+storage_drop_lines + [f'clean()']
    with open(input_file_name,'w') as wr:
        wr.write('\n'.join(content))
    return input_file_name, storage_file_name




def execute_psi4_shell(input_file, psi4_executable=None ):
    
    import subprocess
    if psi4_executable is None:
        psi4_executable=os.path.join( os.path.dirname(os.path.realpath(sys.executable)),  'psi4' )
        if not os.path.isfile(psi4_executable): raise Exception(f"Cannot find 'psi4' executable in same folder as  {sys.executable}")
    #else:
    #    assert os.path.isfile(conda_exc), f"The conda executable does not exist: {conda_exc}"
    #    psi4_executable=f"{conda_exc} run -n {conda_env} psi4"

    output_file=input_file.split('.')[0]+'.psi4out'
    cmd=f"{psi4_executable} {input_file} {output_file}"
    print(f"Excuting \'{cmd}\' in shell")
    start=time.time()
    psi4_run=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = psi4_run.communicate()
    return_code=psi4_run.returncode
    if return_code != 0:
        quit(f"Psi4 finished with return code {return_code}:\n   stdout={stdout.decode()}  error={stderr.decode()}")
    else:
        print(f"Psi4 finished successfully (after {time.time()-start:.3f} seconds)")
    
    return {
        'psi4_output_file': output_file,
        'stdout': stdout.decode(),
        'stderr': stderr.decode(),
        'return_code': return_code,
    }

def psi4_run(
    job,
    func=None,
    basis_set=None,
    xyz_file=None,
    total_charge: int|None=None,
    multiplicity: int|None=None,
    ac_shift: float|None=None,
    #
    execute=False,
    executable=None,
    input_file_name=None,
    do_test=False,
    num_threads=4,
    memory='4GB',
    add_method_tags=True,
):
    if total_charge is None: raise Exception(f"Specify total charge (for safety)")
    if multiplicity is None: raise Exception(f"Specify multiplicity (for safety)")

    try:
        float(memory)
        memory=f"{memory.strip()}GB"
    except:
        assert any([ memory.upper().endswith(x) for x in ['MB','GB'] ]), f"wrong memory {memory}"


    import logging
    logging.basicConfig()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    my_info=logger.info

    
    from .trim_input import psi4_basis_set_escapes, check_basis_known, check_func_known

    if not basis_set is None:
        basis_set=psi4_basis_set_escapes(basis_set)
    
    
    input_and_storage_files=[]
    # psi4 is only needed for execution (if even)
    

    func=check_func_known(func)
    basis_set=check_basis_known(basis_set)

    
    # It test set up auto xyz file
    if not do_test:
        if isinstance(xyz_file, list):
            pass
        elif isinstance(xyz_file, str):
            xyz_file=[xyz_file]
        else:
            raise Exception(f"Did not provide any xyz files!")
    else:
        # Generate dummy xyz
        xyz_file='H2_dummy.xyz'
        dummy_xyz='\n'.join(['2','Dummy file for tests','H 0 0 0','H .74 0 0'])
        with open(xyz_file,'w') as wr:
            wr.write(dummy_xyz)
        xyz_file=[xyz_file]

    for xyz_fi in xyz_file:
        # Get geometry from xyz file
        try:
            atom_types, coordinates = xyz_file_read(xyz_fi)
            charge=total_charge
            multiplicity=multiplicity
            molecule_str='\n'.join([4*' '+x for x in molecule_string(
                atom_types=atom_types, coordinates=coordinates, total_charge=total_charge, multiplicity=multiplicity, units={'LENGTH':'ANGSTROM'}).split('\n') ])
            molecule_name="monomer"
            molecule_lines=[ f"#Define the molecule (definition in molecule string has the advantage of being able to create new molecules reliable)",
                f"monomer_string=\"\"\"\n{molecule_str}\n\"\"\"", f"{molecule_name}=Molecule.from_string(monomer_string)"]
            geom={
                'molecule_name':    molecule_name,
                'charge':           charge,
                'multiplicity':     multiplicity,
                'ac_shift':       ac_shift,
                'atom_types':       atom_types,
                'coordinates':     coordinates,
                'molecule_lines':   molecule_lines,
            }
        except Exception as ex:
            raise Exception(f"Could not process {xyz_fi} as xyz_file: {ex}")

        open_spin=(not multiplicity==1)
        
        # For given geometry get all input files
        for f in func:
            for b in basis_set:
                job_mapper={'polarizabilities':'POLARIZ',
                            'single_point':'SP',
                    }
                assert job in job_mapper.keys(), f"No tag defined for {job}: {job_mapper}"


                jobname=f"{os.path.basename(xyz_fi).split('.')[0]}"
                if add_method_tags:
                    jobname=f"{f}_{b}_{jobname}"

                input_file_name=f"{job_mapper[job].upper()}_{jobname}.psi4inp"
                storage_file_name=f"STORAGE_{job_mapper[job].upper()}_{f}_{b}_{os.path.basename(xyz_fi).split('.')[0]}.yaml"
                input_file_name, storage_file_name = psi4_input_single(
                    job, f, b, geom, open_spin, input_file_name, storage_file_name=storage_file_name, jobname=jobname, 
                    num_threads=num_threads, memory=memory)
                input_and_storage_files.append({'input':input_file_name, 'storage':storage_file_name})
                
        for di in input_and_storage_files:
            input_file_name=di['input']
            storage_file=di['storage']
            if execute:
                execute_psi4_shell(input_file_name, psi4_executable=executable)
            else:
                print(f"Wrote the {job} instructions to {input_file_name}")
        return input_and_storage_files