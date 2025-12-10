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
    POLARIZ = 'POLARIZ'
    SP = 'SP'
    POLARIZ_FF = 'POLARIZ_FF'
    OPT = 'OPT'

class psi4_job(myBaseModel):
    jobtag : job_tags
    jobname : str|None=None
    add_method_tags : bool = True
    
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
            polarizabilities=job_tags.POLARIZ,
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
    def __init__(self, **data):
        jobtag_key='jobtag'
        if jobtag_key in data.keys():
            jobtag=data[jobtag_key]
            if jobtag.upper() not in self.__job_mapper__.values():
                if jobtag.lower() in self.__job_mapper__.keys():
                    data[jobtag_key]=self.__job_mapper__[jobtag.lower()]
                else:
                    raise Exception(f"Job name {jobtag} not recognized! Choose from {list(self.__job_mapper__.keys())} or {list(self.__job_mapper__.values())}")
        super().__init__(**data)
        if self.jobname is None:
            self.jobname=self.make_job_name()
        if self.files.storage is None:
            self.files.storage=self.make_storage_name()
        if self.files.input is None:
            self.files.input=f"{self.jobname}.psi4inp"



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

        return lines