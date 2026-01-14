from . import *
from .geom import geometry, geometry_wrapper
from .machine_settings import machine_settings, my_machine_settings
from .method import psi4_method


class variable_names(BaseModel):
    storage: str =f"storage" # internally used variable name
    molecule_name: str = f"the_molecule" # call the two make clear that this is an instance
    molecule_string: str = f"monomer_string"
    storage_file_name: str = f"storage_file" # internally used variable name
class file_names(BaseModel):
    storage: str | None=None # The file that will be generated after the run
    input: str | None=None # The input file name
    fchk: str|None=None

class job_tags(str, Enum):
    POLARIZ_LR = 'POLARIZ_LR'
    SP = 'SP'
    POLARIZ_FF = 'POLARIZ_FF'
    OPT = 'OPT'


example_geometry_entry=dict(
                input='<provide path to xyz files here>',
                charge=0,
                multiplicity=1,
                ac_shift=None,
)

from .job_settings import (
    pol_ff_options,
    optim_options
)

class psi4_job(controller_model):
    __example_values__=dict(
        geom=example_geometry_entry,
        method=psi4_method.example(),
    )
    jobtag : job_tags
    jobname : str|None=None
    add_method_tags : bool = True
    job_settings: pol_ff_options| optim_options | None =None 
    
    geom : geometry_wrapper
    method: psi4_method 

    execute : bool = False
    executable : str|None =None
    machine_settings : my_machine_settings = my_machine_settings()

    var_names : variable_names = variable_names()
    files : file_names = file_names()

    do_moments : bool = True

    _available_job_tags= job_tags
    __job_mapper__=dict(
            polarizability_linear_response=job_tags.POLARIZ_LR,
            single_point=job_tags.SP,
            polarizability_finite_field=job_tags.POLARIZ_FF,
            optimization=job_tags.OPT,
    )
    @classmethod
    def escape(cls,x):
        return x
        #return x.replace('(','\\(').replace(')','\\)')
    def make_job_name(self):
        jobname=f"{os.path.basename(self.geom.ifi).split('.')[0]}"
        if self.add_method_tags:
            jobname=f"{self.jobtag}_{self.escape(self.method.method_tag)}_{self.escape(self.method.basis_set)}_{jobname}"
        return jobname
    def make_storage_name(self):
        return f"STORAGE_{self.jobname}.yaml"
    
    @classmethod
    def get_settings_for_job_type(cls, job_type):
        mapper={
            job_tags.POLARIZ_FF: pol_ff_options,
            job_tags.OPT: optim_options,
        }
        if not job_type in mapper:
            raise NotImplementedError(f"Job options for job type {job_type} not implemented yet, please provide them manually!")
        return mapper[job_type]()
    @classmethod
    def recognized_job_types(cls, jobtag):
        if jobtag.upper() not in cls.__job_mapper__.values():
            if jobtag.lower() in cls.__job_mapper__.keys():
                return cls.__job_mapper__[jobtag.lower()]
            else:
                raise Exception(f"Job name {jobtag} not recognized! Choose from {list(self.__job_mapper__.keys())} or {list(self.__job_mapper__.values())}")
    def __init__(self, **data):
        jobtag_key='jobtag'
        if jobtag_key in data.keys():
            data[jobtag_key]=self.recognized_job_types( data[jobtag_key] )
        super().__init__(**data)
        if not self.job_settings:
            self.job_settings=self.get_settings_for_job_type(self.jobtag)

            # Use the implicit one
        if isinstance(self.job_settings, pol_ff_options):
            can_do_moments=self.job_settings.moments_possible(self.method.scf_type, self.method.method)
            if self.job_settings.do_moments and not can_do_moments: # if moments requested but not possible
                warn(f"Requested calculation of molecular moments but not possible with method {self.method.method} and scf_type {self.method.scf_type}, setting do_moments to {can_do_moments}")
                self.job_settings.do_moments=can_do_moments

        if self.jobname is None:
            self.jobname=self.make_job_name()
        if self.files.storage is None:
            self.files.storage=self.make_storage_name()
        if self.files.input is None:
            self.files.input=f"{self.jobname}.psi4inp"
    @classmethod
    def example(cls,**kwargs):
        dat=super().example(**kwargs)

        dat['job_settings']=cls.get_settings_for_job_type( cls.recognized_job_types(dat['jobtag']) ).example()  

        # Remove run data
        dat.pop('var_names',None)
        dat.pop('files',None)

        return dat



allowed_job_choices=list(psi4_job.__job_mapper__.keys())

@val_call
def store(job:psi4_job,path:tuple, target_var:str, update:bool=False):
    storage_var=job.var_names.storage
    the_path=''.join([f"[\'{p}\']" for p in path])
    return [ f"{storage_var}{the_path}"+
        ( f"= {target_var}" if not update else
            f".update({target_var})"
        ) ]

class psi4_job(psi4_job):
    def method_lines(self)->List[str]:
        lines=self.method.method_definitions( storage_var=self.var_names.storage )

        method=self.method
        if method.do_grac:
            lines+=method.make_grac_lines(molecule_name=self.var_names.molecule_name, molecule_string=self.var_names.molecule_string)

        if method.do_ip_tuning:
            lines+= method.get_ip_tuned_omega(storage_var=self.var_names.storage)

        new_lines,opt_var=self.method.method_setting_lines(self.geom, return_opts_var=True)
        lines+=new_lines+[
            f"# Store method options",
            *store( self, ('settings',), opt_var, update=True )
        ]
        from qcp_bsx_interface.augmentation.produce_augmentation import requires_augmentation, produce_augmentation, make_psi4_block
        if requires_augmentation(self.method.basis_set):
            lines+= [ 
                f""
                f"# Basis set {self.method.basis_set} is not standard in psi4 distribution but requires augmentation through basis_set_exchange api" ]+[
                *make_psi4_block(produce_augmentation(self.method.basis_set, elements=self.geom.atom_types), basis_name=self.method.basis_set)
            ]

        return lines