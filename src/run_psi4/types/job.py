from . import *
from .geom import geometry, geometry_wrapper
from .machine_settings import machine_settings, my_machine_settings
from .method import psi4_method

### Make variable file registry
# class variables(BaseModel):
#     the_vars: Dict[str,str]



class variable_names(BaseModel):
    storage: str =f"storage" # internally used variable name
    molecule_name: str = f"the_molecule" # call the two make clear that this is an instance
    molecule_string: str = f"monomer_string"
    storage_file_name: str = f"storage_file" # internally used variable name
class file_names(BaseModel):
    storage: str | None=None # The file that will be generated after the run
    input: str | None=None # The input file name
    fchk: str|None=None

class psi4_job(myBaseModel):
    jobtag : str
    jobname : str|None=None
    add_method_tags : bool = True
    
    geom : geometry_wrapper
    method_tag : str
    basis_set : str
    method: psi4_method | None=None

    execute : bool = False
    executable : str|None =None
    machine_settings : my_machine_settings = my_machine_settings()

    var_names : variable_names = variable_names()
    files : file_names = file_names()
    
    __job_mapper__=dict(
            polarizabilities='POLARIZ',
            single_point='SP',
            polarizability_finite_field='POLARIZ_FF',
            optimization='OPT',
    )
    @classmethod
    def escape(cls,x):
        return x
        #return x.replace('(','\\(').replace(')','\\)')
    def make_job_name(self):
        jobname=f"{os.path.basename(self.geom.ifi).split('.')[0]}"
        if self.add_method_tags:
            jobname=f"{self.jobtag}_{self.escape(self.method_tag)}_{self.escape(self.basis_set)}_{jobname}"
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

        if self.method is None:
            self.method=psi4_method(self.method_tag
            )


allowed_job_choices=list(psi4_job.__job_mapper__.keys())
