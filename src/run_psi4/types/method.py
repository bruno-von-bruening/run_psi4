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
    ret['method_corrections']=[]
    for comp in components:
        if check_grac(comp):
            ret['do_grac']=True
        elif check_ip_tuning(comp):
            ret['do_ip_tuning']=True
        elif comp.lower() in ['v','d3bj','d','lrc','a','l']:
            ret['method_corrections']+=[comp.lower()]
        else:
            method_tags.append(comp)
    
    # Check the remainder (should be one component that is the method string)
    if len(method_tags)==0:
        raise Exception(f"Cannot process method string {method}, no valid method found")
    elif len(method_tags)==1:
        ret['method']='-'.join([method_tags[0]]+ret.pop('method_corrections'))
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

def method_setup_lines(method, geom, return_opts_var=False):

    meth_fl=p4ut.method_types
    meth=p4ut.get_method_flags(method.method)



    the_options=dict(
        basis=method.basis_set_var,
        e_convergence=(1.e-6, f"psi4.default"),
        d_convergence=(1.e-6, f"psi4.default"),
        freeze_core=(method.freeze_core, f"Freeze core is significantly cheaper and little impact on e.g. polarizabilities"),
    )

    if meth==p4ut.method_types.MP:      
        if method.method.lower()=='mp2':
            func=p4ut.make_mp2_options
        elif method.method.lower()=='mp3':
            func=p4ut.make_mp3_options
        elif method.method.lower()=='mp4':
            func=p4ut.make_mp4_options
        else: raise Exception(f"Unkown method {method.method}")
    elif meth==p4ut.method_types.CC:    func=p4ut.make_cc_options
    elif meth==p4ut.method_types.CI:    func=p4ut.make_ci_options
    elif meth==p4ut.method_types.HF:    func=p4ut.make_hf_options
    elif meth==p4ut.method_types.DFT:   func=p4ut.make_dft_options
    else:
        raise ValueError(f"Method {method.method} not recognized!")

    if method.reference.lower()=='auto':
        reference=(
                '\'uhf\'' if (geom.open_shell or method.unrestricted) and not method.do_grac  # uhf not possible with grac!
                else '\'rhf\''
        )
    else:
        reference=f"'{method.reference}'"
    the_options.update( func(scf_type=method.scf_type), reference=reference )

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
    def form(v):
        def main_form(v):
            return v.strip('{}') if isinstance(v,str) else str(v)
        if isinstance(v,tuple):
            return ( main_form(v[0]), ' # '+v[1] )
        else:
            return ( main_form(v), '' )


    options_set=[ 
        "{}{:<15} = {:<15},{}".format(
            indent, k.strip(), *form(v) 
        ) 
        for k,v in the_options.items() 
    ]
    # also print the options into the psi4 output 
    # options_printout=[
    #     f"print_out(\"{indent}{k}={{}}\\n\".format({v}))" for k,v in the_options.items()
    # ]
    info_lines=[f"# Information about the settings printed to output file", 
                f"print_out(f\"INFORMATION about SET METHOD\\n\")",
                f"print_out('\\n'.join([f\"{indent}{{k}}={{v}}\" for k,v in opts.items()])) ",
    ]
    #info_lines+=options_printout

    # Create the whole section
    method_setup_lines=[
        '# Setting the options for the calculation',
        "opts=dict("] + options_set + [")",
        'clean_options() # Set options to default (avoid bad surprises)',
        "set_options(opts)",
    ]

    lines= method_setup_lines+info_lines
    if return_opts_var:
        return lines, 'opts'
    else:
        return lines

class psi4_method(myBaseModel):
    #__example_values__=dict(
    method_tag:str
    method: str
    basis_set: str
    freeze_core: bool = True
    scf_type : str='default'
    reference: Literal['AUTO','RHF','UHF','ROHF']='AUTO'

    # scf_type: str='df'
    do_ip_tuning: bool = False
    do_grac: bool = False

    # @classmethod
    # def get_default_scf_type(cls, method):

    @classmethod
    def example(cls):
        # Check that enough mandatory information is provided!
        dummy_info=dict(
            method_tag='HF',
            basis_set='sto-3g',
            freeze_core=True,
            scf_type='default',
        )
        try:
            cls(**dummy_info)
        except Exception as e:
            raise RuntimeError(f"Cannot create example psi4_method: {e}")

        dummy_info['method_tag']='<method[-grac,-iptuned]>'
        return dummy_info

    def __init__(self, *args, **kwargs):
        if len(args)==0:
            if 'method_tag' in kwargs.keys():
                new_kwargs=process_method(kwargs['method_tag'])
            else:
                new_kwargs={}
        elif len(args)==1:
            new_kwargs=process_method(args[0])
            if 'method_tag' in kwargs.keys():
                assert kwargs['method_tag']==args[0], f"Conflict in passed method_tag: {kwargs['method_tag']} vs {args[0]}"
            else:
                kwargs['method_tag']=args[0]
        else:
            raise Exception(f"Passed more than one positional argument to psi4_method: {args}")
        for k,v in new_kwargs.items():
            if k in kwargs:
                if kwargs[k]!=v:
                    warn(f"Conflict in passed arguments to psi4_method: {k} has value {kwargs[k]} from kwargs and {v} from positional argument, will take one from string")
                else: pass
            else:
                kwargs[k]=v
        super().__init__(**kwargs)
class psi4_method(psi4_method):
    func_var: str='method'
    basis_set_var: str='basis_set'
    grac_shift_var: str="grac_shift"
    omega_var: str='omega'
    freeze_core : bool = True

    unrestricted : bool = False

    @property
    def can_do_gradient(self)->bool:
        forbidden=[ r"cc.*\(.*\)", r"mp[2-9]", r"wpbeh-lrc",r"b97" ]
        return False if any(
            re.match(x, self.method.lower()) and self.scf_type.lower()=='conv'
            for x in forbidden
        ) else True
    @val_call
    def method_definitions(self, storage_var:str) -> List[str]:
        assert self.basis_set!='', "Basis set cannot be empty string!"
        return [ x.strip() for x in 
            f"""\
            # Description of the method (indicate correction in comments)
            {self.func_var}=\"{self.method}\"
            {self.basis_set_var}='{self.basis_set}'
            # IP_TUNING: {self.do_ip_tuning}
            # GRAC correction: {self.do_grac}
            {storage_var}['settings'].update({{\'method\': method, \'grac_correction\':{self.do_grac}, \'ip_tuning\':{self.do_ip_tuning}}})"""
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
    def method_setting_lines(self, geom: geometry, return_opts_var:bool=False):
        return method_setup_lines(self, geom, return_opts_var=return_opts_var)






