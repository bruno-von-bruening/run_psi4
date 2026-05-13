from . import *
from run_psi4.types.geom import geometry_model
from run_psi4.config.auto_options import method_types, get_method_flags, make_cc_options,make_ci_options,make_dft_options,make_hf_options,make_mp_options, make_mp_options_gen
from qcp_bsx_interface.bs_modification.process_basis import bs_requires_modification
from qcp_bsx_interface.bs_modification.process_basis import modify_bs,make_psi4_block

def check_grac(comp):
    return  comp.upper() in ['GRAC']
def check_ip_tuning(comp):
    ip_tags=['IPTUNED','IPTUNING']
    return comp.upper() in ip_tags
def decorate_string(scf_type:str):
    return "'{}'".format(scf_type.strip('\'\"'))

@val_call
def process_method(method:str) -> Dict[str,str|bool]:
    # Check for grac
    separator='-'
    components=method.upper().split(separator)

    ret={}

    method_tags=[]
    ret['method_corrections']=[]
    ret['method-prefix']=[]
    
    for comp in components:
        if check_grac(comp):
            ret['do_grac']=True
        elif check_ip_tuning(comp):
            ret['do_ip_tuning']=True
        elif ( comp.lower() in ['cam','lrc','core','dsd','hjs','t','tuned']):
            ret['method-prefix']+=[comp.lower()]
        elif ( comp.lower() in ['v','d','lrc','a','l','nl','gga1','k','vwn','2x','hx','so','rc04','sx','x','sol','dh','fc','fco','fcfo','fl'] or 
            comp.lower().startswith('d3')     or
            bool(re.match(r"[0-9]+[a-z]?", comp.lower()))
        ):
            ret['method_corrections']+=[comp.lower()]
        else:
            method_tags.append(comp)
    
    # Check the remainder (should be one component that is the method string)
    if len(method_tags)==0:
        raise Exception(f"Cannot process method string {method}, no valid method found")
    elif len(method_tags)==1:
        ret['method']='-'.join( ret.pop('method-prefix') +[method_tags[0]]+ret.pop('method_corrections'))
    else:
        raise Exception(f"Cannot process method string {method}, multiple method tags found: {method_tags}")
    return ret



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

    meth=method.method_type

    # Basis options
    the_options=dict(
        basis=method.basis_set_var,
        freeze_core=(method.freeze_core, f"Freeze core is significantly cheaper and little impact on e.g. polarizabilities"),
        nat_orbs=(method.fno, f"Freeze natural orbitals"),
        occ_tolerance=(method.fno_thres if method.fno else 0, f"Freeze natural orbitals threshold"),
        cholesky_tolerance=(method.cd_thres, f"Cholesky decomposition threshold (only relevant for type='cd')"),
    )
    # On frozen natural orbitals: https://psicode.org/psi4manual/1.9.x/fnocc.html

    dft_options=None # the default # printed separately
    if meth==method_types.MP:      
        func=make_mp_options_gen(method.method)
    elif meth==method_types.CC:    func=make_cc_options
    elif meth==method_types.CI:    func=make_ci_options
    elif meth==method_types.HF:    func=make_hf_options
    elif meth==method_types.DFT:   
        func=make_dft_options
        dft_options=method.dft_opts.form_dft_options()
    else:
        raise ValueError(f"Method {method.method} not recognized!")

    conv_opts=method.convergence.model_dump()
    if meth in [ method_types.HF, method_types.DFT ]:
        conv_opts.pop('r_convergence', None) # Remove r_convergence for HF and DFT since it is not relevant for these methods

    if method.reference.lower()=='auto':
        if method.dft_opts.do_grac and geom.open_shell:
            raise Exception(f"Cannot do GRAC correction with open shell systems and automatic reference detection (UHF not possible with GRAC), please set reference explicitly!")
        reference=(
                '\'uhf\'' if (geom.open_shell) else '\'rhf\''
        )
    else:
        reference=f"'{method.reference}'"
    the_options.update( func(scf_type=method.scf_type), reference=(reference.upper() if not reference is None else None) )

    df_options={}
    for ty in ['cc','mp2','scf']:
        if f"{ty}_type" in the_options.keys():
            if the_options[f"{ty}_type"].strip('\'\"').lower() in ['df']:
                df_options[f"df_basis_{ty}"]=getattr(method.df_opts, f"df_basis_{ty}")


    if method.dft_opts.do_ip_tuning:
        dft_options.update(dict(
            DFT_omega=f"{{{method.omega_var}}}"
        ))
    if method.dft_opts.do_grac:
        dft_options.update(dict(
            DFT_GRAC_SHIFT=var_str(f"{method._grac_shift_var}"),
        ))

    # Print the options into the file
    indent=4*' '
    def form(v):
        def main_form(v):
            base=v.strip('{}') if isinstance(v,(str)) else str(v)
            if not issubclass(type(v), var_str) and issubclass(type(v), str):
                base=decorate_string(base)
            return base
        if isinstance(v,tuple):
            return ( main_form(v[0]), ' # '+v[1] )
        else:
            return ( main_form(v), '' )


    def set_options(di):
        return [ 
            "{}{:<15} = {:<15},{}".format(
                indent, k.strip(), *form(v) 
            ) 
            for k,v in di.items() 
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


    # Create the sub options
    set_opts=[]
    for opt, name in zip(
        [ dft_options, df_options, conv_opts],
        [ 'dft','df','conv']
    ):
        if opt:
            set_opts+=[
                f"{name}_opts=dict(", *set_options(opt), ')',
            ]

    df_exclude=[ x for x in df_options.values() if not x is None and bs_requires_modification(x) ]
    df_exclude_str=( f" and not v in {df_exclude} # they are custom so have to be set expliclity" if len(df_exclude)>0 else "" )
    set_opts+=[
        "opts=dict(", 
        *set_options(the_options),
        *([] if not dft_options else [indent+"**dft_opts,"]),
        *([] if not df_options else [indent+"**dict( (k,v) for k,v in df_opts.items() ",
                                     2*indent+f"if not v is None {df_exclude_str}",
                                     indent+f"), "]),
        *([] if not conv_opts else [indent+"**conv_opts,"]),
        ")",
        'clean_options() # Set options to default (avoid bad surprises)',
        "set_options(opts)",
    ]

    method_setup_lines=[
        '# Setting the options for the calculation', *set_opts
    ]

    lines= method_setup_lines+info_lines
    if return_opts_var:
        return lines, 'opts'
    else:
        return lines


class dft_opts_model(controller_model):
    # DFT grid options # https://psicode.org/psi4manual/1.9.x/dft.html
    # one to one psi4 keywords
    radial_scheme: Literal['TREUTLER','BECKE', 'MULTIEXP', 'EM', 'MURA']='TREUTLER' # https://psicode.org/psi4manual/1.9.x/autodoc_glossary_options_c.html#term-DFT_RADIAL_SCHEME-SCF
    radial_points: int=75 # https://psicode.org/psi4manual/1.9.x/autodoc_glossary_options_c.html#term-DFT_RADIAL_POINTS-SCF
    spherical_points: int=302 # # https://psicode.org/psi4manual/1.9.x/autodoc_glossary_options_c.html#term-DFT_SPHERICAL_POINTS-SCF
    spherical_scheme: Literal['LEBEDEV']='LEBEDEV'
    # dft_pruning_scheme: Literal['ROBUST','TREUTLER']|None=None # https://psicode.org/psi4manual/1.9.x/autodoc_glossary_options_c.html#term-DFT_PRUNING_SCHEME-SCF
    # Other options
    do_ip_tuning: bool = False
    do_grac: bool = False
    def form_dft_options(self):
        base_di=self.model_dump(exclude=['do_ip_tuning','do_grac'])
        base_di={ 'dft_'+k:v for k,v in base_di.items() if v is not None }
        return dict(
            **base_di,
            # DFT_pruning_scheme=self.dft_pruning_scheme,
        )
# class psi4_method(myBaseModel):

class df_settings_model(controller_model):
    df_basis_scf: str|None=None
    df_basis_mp2: str|None=None
    df_basis_cc: str|None=None
    def __init__(self, input:None|str=None, **kwargs):
        if input is not None:
            for k in ['df_basis_scf','df_basis_mp2','df_basis_cc']:
                if not k in kwargs.keys():
                    kwargs[k]=input
        super().__init__(
            **kwargs
        )

class convergence_model(controller_model):
    e_convergence:  float=1.e-6 # Energy convergence
    d_convergence:  float=1.e-6 # Density matrix convergence
    r_convergence:  float=1.e-6 # RMS convergence of amplitudes
    maxiter:        int=100 # Default psi4 value

ReferenceType = Annotated[
    Literal['AUTO', 'RHF', 'UHF', 'ROHF'],
    BeforeValidator(lambda v: v.upper() if isinstance(v, str) else v)
]
# Make method uppercase!
MethodType=Annotated[ str, BeforeValidator(lambda v: v.upper() if isinstance(v, str) else v) ]

class psi4_method(controller_model):
    #__example_values__=dict(
    __method_types__=method_types
    method_tag: MethodType
    method: MethodType
    basis_set: str
    freeze_core: bool = False
    fno: bool = False # Freeze natural orbitals threshold (any value not None or zero activates it)
    fno_thres: float=1.e-6
    # The 

    scf_type : str='default'
    cd_thres : float=1.e-4 
    reference: ReferenceType = 'AUTO'

    dft_opts: dft_opts_model=dft_opts_model()
    df_opts: df_settings_model=df_settings_model()
    convergence: convergence_model=convergence_model()

    # scf_type: str='df'

    # @classmethod
    # def get_default_scf_type(cls, method):
    __example_values__=dict(
        method_tag='<method[-grac,-iptuned]>',
        basis_set='<basis_set>'
    )
    @property
    def method_type(self):
        return get_method_flags(self.method)

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
            raise RuntimeError(f"Cannot create example psi4_method: {e}\n{dummy_info}")

        return super().example(exclude=['method'])

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

        do_grac=new_kwargs.pop('do_grac', False)
        do_ip_tuning=new_kwargs.pop('do_ip_tuning', False)
        for k,v in new_kwargs.items():
            if k in kwargs:
                if kwargs[k]!=v:
                    warn(f"Conflict in passed arguments to psi4_method: {k} has value {kwargs[k]} from kwargs and {v} from positional argument, will take one from string")
                else: pass
            else:
                kwargs[k]=v
        
        scf_type=kwargs.get('scf_type',None)
        if scf_type:
            if scf_type.lower().startswith('df'):
                if not 'df_opts' in kwargs.keys():
                    tag,*args= scf_type.split('_')
                    if len(args)>0:
                        assert tag.lower()=='df', f"scf_type {scf_type} not recognized, expected to start with 'df'"
                        kwargs['scf_type']=tag
                        assert len(args)==1, f"scf_type {scf_type} not recognized, expected format 'df_<type>'"
                        df_sets=df_settings_model(args[0])
                        kwargs['df_opts']=df_sets.model_dump()

        super().__init__(**kwargs)
        if do_grac:
            self.dft_opts.do_grac=True
        if do_ip_tuning:
            self.dft_opts.do_ip_tuning=True
        if self.scf_type.lower()=='default':
            method_type=self.method_type
            if method_type in [method_types.MP, method_types.CC, method_types.CI]:
                self.scf_type='conv'
            elif method_type in [method_types.HF, method_types.DFT]:
                self.scf_type='pk'
            else: 
                raise ValueError(f"Method {self.method} not recognized!")
        elif self.scf_type.lower()=='conv' and self.method_type in [method_types.HF, method_types.DFT]:
            warn(f"SCF type 'conv' is not a valid key for HF and DFT methods! Automatically changing to 'pk'.")
            self.scf_type='pk' # was important when inherting settings from e.g. ccsd
    def necessary_df_bs(self)->List[str]:
        if self.method_type in [method_types.HF, method_types.DFT]:
            req=['scf']
        elif self.method_type==method_types.MP:
            req=['scf',self.method_tag.lower()]
        elif self.method_type==method_types.CC:
            req=['scf','mp2','cc']
        else:
            raise NotImplementedError(f"Method type {self.method_type} not recognized for determining necessary df basis sets!")
        keys=set( f"df_basis_{k}" for k in req )
        return dict( (k,getattr(self.df_opts, k)) for k in keys )

    
    def basis_set_lines(self, atom_types: List[str])->List[str]:
        # if requires_augmentation(self.method.basis_set):
        if bs_requires_modification(self.basis_set):
            basis_dict=modify_bs(self.basis_set, elements=atom_types)
            bas_lines= [ 
                f""
                f"# Basis set {self.basis_set} is not standard in psi4 distribution but requires augmentation through basis_set_exchange api" ]+[
                *make_psi4_block(basis_dict, basis_name=f"{self.basis_set}")
            ]
        else:
            bas_lines=[]
        
        df_bs_lines=[]
        if self.scf_type.lower()=='df':
            df_basis_sets=self.necessary_df_bs()
            df_bs_needed=set()
            for df_type,bs_nam in df_basis_sets.items():
                if not bs_nam is None and bs_requires_modification(bs_nam):
                    df_bs_needed.add(bs_nam)
            if len(df_bs_needed)>0:
                # psi4 requires expliclity section of defining df basis sets
                # for bs_nam in df_bs_needed:
                #     basis_dict=modify_bs(bs_nam, elements=self.geom.atom_types)
                #     df_bs_lines+= [ 
                #         *make_psi4_block(basis_dict, basis_name=f"{bs_nam}")
                #     ]
                # df_bs_lines=[ 
                #     '',
                #     f"# DF basis sets {df_bs_needed} are non standard and will be hence defined specifially below",
                #     *df_bs_lines
                # ]
                for k,v in df_basis_sets.items():
                    if v in df_bs_needed:
                        basis_dict=modify_bs(bs_nam, elements=atom_types)
                        df_bs_lines+= [ 
                            *make_psi4_block(basis_dict, basis_name=f"{k}", is_df_basis=True)
                        ]
        return bas_lines+df_bs_lines

class var_str(str):
    pass
    # def decorate_string(self):
    #     return decorate_string(self)

class psi4_method(psi4_method):
    _func_var: str='method'
    _basis_set_var: var_str='basis_set'
    _grac_shift_var: var_str="grac_shift"
    _omega_var: var_str='omega'
    freeze_core : bool = True
    @property
    def basis_set_var(self):
        return var_str(self._basis_set_var)

    

    # _unrestricted : bool = False

    @property
    def can_do_gradient(self)->bool:
        forbidden=[ r"cc.*\(.*\)", r"mp[2-9]", r"wpbeh-lrc",r"b97" ]
        if re.fullmatch(r"cc.*", self.method.lower()) and self.scf_type.lower() in ['conv','cd']:
            return False
        return False if any(
            re.match(x, self.method.lower()) and self.scf_type.lower() in ['conv','cd']
            for x in forbidden
        ) else True
    @property
    def mos_implicit(self)->bool:
        return self.method_type in [method_types.HF, method_types.DFT]
    @val_call
    def method_definitions(self, storage_var:str) -> List[str]:
        assert self.basis_set!='', "Basis set cannot be empty string!"
        return [ 
            f"# Description of the method (indicate correction in comments)",
            f"{self._func_var}=\"{self.method}\"",
            f"{self.basis_set_var}='{self.basis_set}'",
            f"{storage_var}['settings'].update({{\'method\': method, \'grac_correction\':{self.dft_opts.do_grac}, \'ip_tuning\':{self.dft_opts.do_ip_tuning}}})"
        ]

    @val_call
    def compute_grac_shift(self,
        func_var, basis_set_var, molecule_name, molecule_string, charge=0, spin=1
    ) ->  List[str] :# Tuple[List[str],str]:
        from run_psi4.main.input_writer import get_gs_energy, get_ionization_energy
        neut_string, E_neut, wfn_neut = get_gs_energy(func_var, basis_set_var, molecule_name=molecule_name)
        ion_string, E_ion, wfn_ion    = get_ionization_energy(func_var, basis_set_var, molecule_name=molecule_name, molecule_string=molecule_string)

        # Compute the difference -> grac_shift
        string=[ 
                '# computing grac shift as IP+HOMO',
                f"E_homo = float({wfn_neut}.epsilon_a_subset('AO', 'ALL').np[{wfn_neut}.nalpha()-1])",
                f"IP = float({E_ion}) - float({E_neut})",
                f"{self._grac_shift_var} = IP + E_homo",
        ]
        store_string=[
            f"# Store the values"
            f"storage['results'].update({{ 'HOMO': E_homo, 'IP': IP , 'AC_SHIFT': {self._grac_shift_var} }})",
        ]

        grac_string=neut_string+ion_string+string+store_string
        return grac_string
    
    @val_call
    def provide_grac(self, ac_shift:float):
        assert self.dft_opts.do_grac, "Requested to make grac string but was internally set false (something is wrong in code logic)"
        grac_string=[f"# Using user provided AC shift of {ac_shift} (assumed to be in Hartree)"]
        grac_lines, grac_var=set_grac_shift(ac_shift)
        grac_string+=grac_lines
        self.grac_shift_var=grac_var
        return 
    @val_call
    def make_grac_lines(self,  molecule_name:str|None=None, molecule_string:str|None=None) -> List[str] :
        assert self.dft_opts.do_grac, "Requested to make grac string but was internally set false (something is wrong in code logic)"
        grac_string=[f"# No grac shift provided, computing HOMO and IP to get AC shift"]
        assert molecule_name is not None, f"Provide moleucle name for grac shift calculation!"
        assert molecule_string is not None, f"Provide moleucle string for grac shift calculation!"
        grac_string+=self.compute_grac_shift(self._func_var,self.basis_set_var, molecule_name, molecule_string)
        return grac_string
    @val_call
    def method_setting_lines(self, geom: geometry_model, return_opts_var:bool=False):
        return method_setup_lines(self, geom, return_opts_var=return_opts_var)






