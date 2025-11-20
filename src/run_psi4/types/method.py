from . import *
from run_psi4.types.geom import geometry

def check_grac(comp):
    return  comp.upper() in ['GRAC']
def check_ip_tuning(comp):
    ip_tags=['IPTUNED','IPTUNING']
    return comp.upper() in ip_tags

@val_call
def process_method(method:str) -> Dict[str,str|bool]:
    # Check for grac
    separator='-'
    components=method.upper().split(separator)

    ret={}

    method_tags=[]
    for comp in components:
        if check_grac(comp):
            ret['do_grac']=True
        elif check_ip_tuning(comp):
            ret['do_ip_tuning']=True
        else:
            method_tags.append(comp)
    
    # Check the remainder (should be one component that is the method string)
    if len(method_tags)==0:
        raise Exception(f"Cannot process method string {method}, no valid method found")
    elif len(method_tags)==1:
        ret['method']=method_tags[0]
    else:
        raise Exception(f"Cannot process method string {method}, multiple method tags found: {method_tags}")
    return ret


@val_call
def compute_grac_shift(
    func_var, basis_set_var, molecule_name, molecule_string, charge=0, spin=1
) -> Tuple[List[str],str]:
    from run_psi4.main.input_writer import get_gs_energy, get_ionization_energy
    neut_string, E_neut, wfn_neut = get_gs_energy(func_var, basis_set_var, molecule_name=molecule_name)
    ion_string, E_ion, wfn_ion    = get_ionization_energy(func_var, basis_set_var, molecule_name=molecule_name, molecule_string=molecule_string)

    # Compute the difference -> grac_shift
    ac_var=f"AC_shift"
    string=[ 
            '# computing grac shift as IP+HOMO',
            f"E_homo = float({wfn_neut}.epsilon_a_subset('AO', 'ALL').np[{wfn_neut}.nalpha()-1])",
            f"IP = float({E_ion}) - float({E_neut})",
            f"{ac_var} = IP + E_homo",
    ]
    store_string=[
        f"# Store the values"
        f"storage['results'].update({{ 'HOMO': E_homo, 'IP': IP , 'AC_SHIFT': {ac_var} }})",
    ]

    grac_string=neut_string+ion_string+string+store_string
    return grac_string,ac_var

@val_call
def set_grac_shift(
    ac_var:str|float
) -> Tuple[List[str],str]:
    grac_var='the_grac_shift'
    return [x.strip() for x in f"""
        set_options({{'dft_grac_shift': {ac_var}}})
        {grac_var}=get_option('SCF','dft_grac_shift')
        print_out(f\"GRAC SHIFT SET TO: {{{grac_var}}}\")
        assert {grac_var}>0.0, f\"The grac shift is not larger zero!\"
        """.split('\n')
    ],grac_var

def method_setup_lines(method, geom):

    meth_flags=p4ut.get_method_flags(method.method)


    the_options=dict(
        basis=method.basis_set_var,
    )

    # lines+=[f"psi4.set_options({{",
    #  f"\"scf_type\": \"df\",",
    # f"\"e_convergence\": 1e-12,",      # energy convergence
    # f"\"d_convergence\": 1e-12,",       # density matrix convergence
    # f"}})"]

    # MP2 related keywords
    if meth_flags['mp2']:
        the_options.update( p4ut.make_mp2_options() )
    elif meth_flags['cc']:
        the_options.update( p4ut.make_cc_options() )
    elif meth_flags['hf']:
        the_options.update( p4ut.make_hf_options(open_spin=geom.open_shell) )
    else:
        the_options.update( p4ut.make_dft_options(open_spin=geom.open_shell) )

    if method.do_ip_tuning:
        the_options.update(dict(
            DFT_omega=f"{{{method.omega_var}}}"
        ))
    if method.do_grac:
        the_options.update(dict(
            DFT_GRAC_SHIFT=f"{method.grac_shift_var}",
        ))

    # Print the options into the file
    indent=4*' '
    options_set=[ 
        "{indent}{k:<15} : {v:<15} ,".format(
            indent=indent, k='\''+k+'\'', v=v.strip('{}')
        ) 
        for k,v in the_options.items() 
    ]
    # also print the options into the psi4 output 
    options_printout=[
        f"print_out(\"{indent}{k}={{}}\\n\".format({v}))" for k,v in the_options.items()
    ]
    info_lines=[f"# Information about the settings printed to output file", f"print_out(f\"INFORMATION about SET METHOD\\n\")"]
    info_lines+=options_printout

    # Create the whole section
    method_setup_lines=[
        '# Setting the options for the calculation',
        'clean_options() # Set options to default (avoid bad surprises)',
        "set_options({"] + options_set + ["})"]
    
    return method_setup_lines+info_lines

class psi4_method(myBaseModel):
    method: str
    do_ip_tuning: bool = False
    do_grac: bool = False
    def __init__(self, *args, **kwargs):
        if len(args)==0:
            pass
        elif len(args)==1:
            new_kwargs=process_method(args[0])
            for k,v in new_kwargs.items():
                if k in kwargs:
                    if kwargs[k]!=v:
                        warn(f"Conflict in passed arguments to psi4_method: {k} has value {kwargs[k]} from kwargs and {v} from positional argument, will take one from string")
                    else: pass
                else:
                    kwargs[k]=v
        else:
            raise Exception(f"Passed more than one positional argument to psi4_method: {args}")
        super().__init__(**kwargs)
class psi4_method(psi4_method):
    func_var: str='method'
    basis_set_var: str='basis_set'
    grac_shift_var: str="grac_shift"
    omega_var: str='omega'
    @val_call
    def method_definitions(self, storage_var:str, basis_set:str) -> List[str]:
        assert basis_set!='', "Basis set cannot be empty string!"
        return [ x.strip() for x in 
            f"""\
            # Description of the method (indicate correction in comments)
            {self.func_var}=\"{self.method}\"
            {self.basis_set_var}='{basis_set}'
            # IP_TUNING: {self.do_ip_tuning}
            # GRAC correction: {self.do_grac}
            {storage_var}['settings'].update({{\'method\': method, 'basis_set': basis_set, \'grac_correction\':{self.do_grac}, \'ip_tuning\':{self.do_ip_tuning}}})"""
            .split('\n')
        ]

    @val_call
    def make_grac_lines(self, ac_shift:float|None=None, molecule_name:str|None=None, molecule_string:str|None=None) -> List[str] :
        assert self.do_grac, "Requested to make grac string but was internally set false (something is wrong in code logic)"
        if ac_shift is not None:
            grac_string=[f"# Using user provided AC shift of {ac_shift} (assumed to be in Hartree)"]
            grac_lines, grac_var=set_grac_shift(ac_shift)
            grac_string+=grac_lines
            self.grac_shift_var=grac_var
        else:
            grac_string=[f"# No grac shift provided, computing HOMO and IP to get AC shift"]
            assert molecule_name is not None, f"Provide moleucle name for grac shift calculation!"
            assert molecule_string is not None, f"Provide moleucle string for grac shift calculation!"
            grac_lines, grac_var=compute_grac_shift(self.func_var,self.basis_set_var, molecule_name, molecule_string)
            self.grac_shift_var=grac_var
            grac_string+=grac_lines
        return grac_string
    @val_call
    def method_setting_lines(self, geom: geometry):
        return method_setup_lines(self, geom)






