from . import *
from run_psi4.types.geom import geometry, geometry_wrapper
from run_psi4.types.machine_settings import machine_settings, my_machine_settings
from run_psi4.types.job import allowed_job_choices, psi4_job

from run_psi4.execution_and_recovery.execute import execute_psi4_shell
from run_psi4.input_file_generation.main import  psi4_input_single

    # total_charge    : int|None=None,
    # multiplicity    : int|None=None,
    # ac_shift        : float|None=None,


from .trim_input import psi4_basis_set_escapes, check_basis_known, check_func_known
def check_basis_set(v):
    if isinstance(v, str):
        try:
            basis_set=psi4_basis_set_escapes(v)
            basis_set=check_basis_known(basis_set)
            return basis_set
        except Exception as ex:
            raise ValueError(f"Basis set {v} not recognized!: {ex}")
    else: raise ValueError(f"Basis set has to be a string!")

my_method=Annotated[str, BeforeValidator(check_func_known)]
my_basis_set=Annotated[str, BeforeValidator(check_basis_set)]

@val_call
def psi4_run(
    jobtag             : str,
    geom            : List[geometry_wrapper]|geometry_wrapper|None=None,
    method          : List[my_method]|my_method|None=None,
    basis_set       : List[my_basis_set]|my_basis_set|None=None,
    #
    execute=False,
    executable=None,
    input_file_name=None,
    do_test=False,
    machine_settings : my_machine_settings = machine_settings(),
    add_method_tags : bool = True,
    do_mom: bool = True,
    freeze_core: bool = True,
    use_df: bool = True,
    unrestricted: bool = False,
) -> List[Tuple[str,str]]:

    input_and_storage_files=[]

    if not do_test:
        if geom is None:
            raise Exception(f"Did not provide any xyz files!")
    else:
        # Generate dummy xyz
        xyz_file='H2_dummy.xyz'
        dummy_xyz='\n'.join(['2','Dummy file for tests','H 0 0 0','H .74 0 0'])
        with open(xyz_file,'w') as wr:
            wr.write(dummy_xyz)
        geom=geometry(xyz_file, charge=0, multiplicity=1)

    if not isinstance(geom, list):
        geoms = [geom]
    if not isinstance(basis_set, list):
        basis_sets = [basis_set]
    if not isinstance(method, list):
        methods = [method]

    input_and_storage_files=[]
    for geom in geoms:
        for f in methods:
            for b in basis_sets:

                the_job=psi4_job(
                    jobtag=jobtag,
                    geom=geom,
                    method=dict(method_tag=f, freeze_core=freeze_core, use_df=use_df,basis_set=b, unrestricted=unrestricted),
                    execute=execute,
                    executable=executable,
                    machine_settings=machine_settings,
                    do_moments=do_mom,
                )
                the_job = psi4_input_single(
                    the_job
                )
                input_and_storage_files.append( (the_job.files.input, the_job.files.storage) )

                if execute:
                    execute_psi4_shell(the_job.files.input, psi4_executable=executable)
                else:
                    print(f"Wrote the instructions to {the_job.files.input}")

    return input_and_storage_files