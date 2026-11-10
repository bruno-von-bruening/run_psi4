from . import *

from run_psi4.util.import_helper import *
from .action import psi4_action
from .tag_matrix import tag_matrix, tensor_tag_from_indices
from .center_of import set_origin

@val_call
def meassure_multipoles_func(job):
    """ Meassure multipoles for a single point
    This is currently a function to be defined within psi4 script
    tag tells perturbation (e.g. xx for quadrupolar field)

    Ideally this should take a wfn
    """
    lines=[]
    
    wfn_internal='wfn'

    def method_works_with_properties(method):
        return False
        return (not method.lower().startswith('cc')) and (not method.lower()=='ccsd')

    # compute dipole and quadrupole moment
    properties=['DIPOLE','QUADRUPOLE','OCTUPOLE','HEXADECAPOLE']
    properties=[f"MULTIPOLE({i})" for i  in range(1,5)]
    if method_works_with_properties(job.method.method):
        kwargs=dict(
            properties=properties,
            return_wfn=False,
            wfn=wfn_internal,
        )
        lines+=[ psi4_action(job,'properties', **kwargs) ]
    else:
        lines+=[
            f"oeprop = psi4.core.OEProp({wfn_internal})",
            *[ f"oeprop.add('{t}')" for t in properties],
            f"oeprop.compute()",
        ]
    # clean scratch files
    # lines+=[f"psi4.core.clean()"]
        
    def trans_prop(tag):
        mapper={
            'MULTIPOLE(1)': 'DIPOLE',
            'MULTIPOLE(2)': 'QUADRUPOLE',
            'MULTIPOLE(3)': 'OCTUPOLE',
            'MULTIPOLE(4)': 'HEXADECAPOLE',
        }
        return mapper.get(tag, tag)
        if any(job.method.method.lower().startswith(x) for x in ['mp','cc','ci']):
            return f"{job.method.method} {tag}"
        else: 
            return f"SCF {tag}"
    # save the moments under variables describing the moment 
    lines+=[
        # "for name in ['CURRENT QUADRUPOLE', 'MULLIKEN_QUADRUPOLE', 'Q(QM)']:",
        # "    try:",
        # "        print(f\"{name} =\", wfn.variable(name))",
        # "    except KeyError:",
        # "        print(f\"No {name}\")",
        # "all_vars = wfn.variables()",
        # "print('All available wfn variables:')",
        # "for var in all_vars:",
        # "    print(var)",
        f"def my_trans(x):",
        IND+f"fac=1",# if len(x)==4 else 2",
        IND+f"ar=np.array(x).reshape(tuple([3]*int(np.emath.logn(3,x.size))))", # e-
        IND+f"return {tag_matrix.__name__}(ar) ",
        f"new_di= dict()",
        f"[ new_di.update(my_trans(wfn.variable(prop) ) ) for prop in {[ trans_prop(x) for x in properties]} ]",
        f"return new_di",
        # "return dict({ ",
        # *[ind+f"'{trans_prop(prop)}' :  my_trans(wfn.variable('{trans_prop(prop).upper()}')) ," for prop in properties], # the -1 is for the conversion from psi4 convention to usual convention (e.g. dipole is negative since it e-)
        # '})'
    ]

    import inspect
    lines=[ 
        *inspect.getsource(tensor_tag_from_indices).split('\n'),
        *inspect.getsource(tag_matrix).split('\n'),
        f"def eval_moments({wfn_internal}):",
        *[IND+l for l in lines],
    ]
    return lines

class handle_molecular_moments(BaseModel):
    def make_functions(self, job:psi4_job):
        """ Returns the lines to write the funciton that get moments """
        return meassure_multipoles_func(job)
    def set_origin(self,job:psi4_job)-> Tuple[List[str],List[str]]:
        return set_origin(job)
    def eval_moments_to(self, job:psi4_job, the_place):
        return [
            f"# EVAL MOMENTS stored under {job.var_names.main_wfn}",
            f"{the_place}=eval_moments({job.var_names.main_wfn})"
        ] 
    def need_gradient(self,job): 
        moments_for_free        = ( job.method.method_type in [ method_types.HF, method_types.DFT ] )
        if moments_for_free:
            return False
        elif not job.method.can_do_gradient:
            raise ValueError(f"Need gradient to get moments, but {job.method} cannot!")
        else:
            return True
