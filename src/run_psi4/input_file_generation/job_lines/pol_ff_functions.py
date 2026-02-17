from . import *

from run_psi4.types.job import psi4_job, store
from run_psi4.types.method import method_types

# def make_multipolar_field():

# @val_call
# def make_external_field(key, center, d) -> Tuple[str,List[str]]:
def make_external_field(key, center, d):
    """ key: str (e.g. 'x', 'xx' or 'xy' )
    center: np.ndarray center of fields
    d: float separation parameter (how far are the charges from center)
    """

    def get_dir(key):
        mapper={
            'x': np.array([1.,0.,0.]),
            'y': np.array([0.,1.,0.]),
            'z': np.array([0.,0.,1.]),
        }
        return mapper[key]
    
    if len(key)==1:
        vec=get_dir(key)
        arrangement=[
            ( -1, -vec),
            ( 1, vec),
        ]
        q= -d**2/2.
    elif len(key)==2:
        a,b=key[0],key[1]
        vec1, vec2=get_dir(a), get_dir(b)
        q= -d**3/3
        if a==b:
            arrangement=[
                ( 1., -vec1),
                ( 1., vec1),
            ]
        else:
            arrangement=[
                ( +.5, + (vec1 + vec2)),
                ( +.5, - (vec1 + vec2)),
                ( -.5, + (vec1 - vec2)),
                ( -.5, - (vec1 - vec2)),
            ]
    else:
        raise ValueError(f"Unsupported key length: {len(key)}")

    arrangement=[
        (c*q, center+v*d/np.linalg.norm(v)) for c,v in arrangement
    ]


    return dict(
        charges=arrangement,
        center=center,
        separation=d,
    )

import numpy as np

def nuclear_external_point_charge_energy(wfn):
    """
    Nuclear interaction energy with external point charges.

    Parameters
    ----------
    wfn : psi4.core.Wavefunction

    Returns
    -------
    float
        Nuclear–external point charge interaction energy (Hartree)
    """

    # Correction not needed for HF like methods
    method = wfn.name().upper()
    if not any( method.startswith(x) for x in ["MP",'CC']):
        return 0

    # External charges (if none, Psi4 returns empty arrays)
    pot=wfn.external_pot()
    mol = wfn.molecule()
    if pot is None: return 0
    else: return pot.computeNuclearEnergy(mol)


def add_energy_vars_to_loc_store(job,loc_store_name):
    """Return lines that update loc_store (name as string) with correlated energies if available.
    These are additional variables (e.g. CCSD(T) yields also CCSD energy as preliminary)
    """
    out=[]
    if job.method.method_type in [ method_types.CC, method_types.MP ]:
        out += [
            f"try:",
            ind+f"{loc_store_name}.update(E_HF=psi4.variable(\"HF TOTAL ENERGY\"))",
        ]
        if job.method.method_type in [ method_types.CC]:
            out += [
                ind+f"{loc_store_name}.update(E_MP2=psi4.variable(\"MP2 TOTAL ENERGY\"))",
            ]
            if job.method.method.lower() == 'ccsd':
                # nothing extra for plain CCSD here
                pass
            elif job.method.method.lower() == 'ccsd(t)':
                out += [
                    ind+f"{loc_store_name}.update(E_CCSD=psi4.variable(\"CCSD TOTAL ENERGY\"))"
                    # -psi4.variable('CCSD(T) CORRELATIOn ENERGY')) # CCSD is acutally CCSD(T)!", 
                ]
            else:
                raise ValueError(f"Unsupported CC method: {job.method.method}")
        out+=[
            f"except Exception as ex:",
            ind+f"warning=f\"WARNING: Some preliminary energies (e.g. HF or MP2) not available for method {job.method.method}, skipping: {{str(ex)}}\"",
            ind+f"print(warning);print_out(\'CUSTOM \'+warning) ",
        ]
    return out


def get_moment_function(job):
    """ For each perturbation tag:
        1. get the right field
        2. compute the wfn
        3. evaluate moments
        4. store moments in registry
    """
    # first get the wfn
    loc_store='store_di'
    lines=[
        f"{loc_store}=dict()"
    ]

    #mom_registry='mom_registry'
    #pot_registry='pot_registry'
    #pot_tag='pot_tag'
    ext_pot_int='ext_pot'
    delta='delta'

    # 1. get the potential
    lines+=[ #f"{ext_pot_int}={pot_registry}[{pot_tag}]",
            f"{ext_pot_int}=( [ (q*{delta},v) for q,v in {ext_pot_int}['charges'] ]",
            f"    if not ( {ext_pot_int} is None or {ext_pot_int}==[]) else None)",
    ]
    

    # 2. compute wfn
    kwargs=dict(
        external_potentials=ext_pot_int,
        return_wfn=True,
        save_under='E',
    )
    action='energy' # only needed if not scf or KS, otherwise gradient updates density matrix
    lines+=[ 
        f"single_point_time=time.time()",]+[
        *psi4_action(job, action, **kwargs) ]+[
        *store_timing(job, 'single_point', 'single_point_time',add=True),
        f"{loc_store}.update(E={kwargs['save_under']})"
    ]


    # 
    lines += add_energy_vars_to_loc_store(job,loc_store_name=loc_store)


    # lines+=[
    #         f"if not {ext_pot_int} is None:", 
    #         ind+f"E_int=0.0",
    #         ind+f"the_pot=psi4.core.ExternalPotential()",
    #         ind+f"for pot in {ext_pot_int}:",
    #         ind+ind+f"the_pot.addCharge(pot[0], *pot[1])",
    #         ind+f"D = np.asarray(wfn.Da()) + np.asarray(wfn.Db())",
    #         ind+f"Vmat = the_pot.computePotentialMatrix(wfn.basisset())",
    #         ind+f"V = np.asarray(Vmat)",
    #         ind+f"E_int += float(np.einsum(\"ij,ij\", D, V))",
    #         ind+f"for i in range(wfn.molecule().natom()):      # natom() returns number of atoms",
    #         ind+f"    Z = wfn.molecule().Z(i)                  # nuclear charge of atom i",
    #         ind+f"    x, y, z = wfn.molecule().x(i), wfn.molecule().y(i), wfn.molecule().z(i)",
    #         ind+f"    E_int += sum(q / ((x-xq)**2 + (y-yq)**2 + (z-zq)**2)**0.5 for q, (xq, yq, zq) in {ext_pot_int})",
    #         ind+f"print_out(f\"Interaction energy: {{E_int}}\\n\")",
    #         f"else: E_int=0.0",
    # ]

# 4) Contract to get interaction energy

    if any( job.method.method.lower().startswith(x) for x in ['cc','mp']):
        kwargs.update(save_under='grad', ref_wfn='wfn')
        action='gradient' # only needed if not scf or KS, otherwise gradient updates density matrix
        grad_lines=[ 
            f"# For CC and MP2 we need to update the density matrix through a gradient calculation",
            f"time_grad=time.time()",
            f"# we need a gradient calculation for evaluation of density matrix",
            *psi4_action(job, action, **kwargs),
            *store_timing(job, 'grad', 'time_grad', add=True),
            ]
        lines+=[ f"if do_mom:"]+[ ind+l for l in grad_lines ]+[f"else: grad=None"]
    else:
        lines+=[f"grad=None "]
    lines+=[f"{loc_store}.update(grad=grad)"]

    # 3. evaluate moments
    lines+=[ 
            '',
            #f"mom= {{'DIPOLE':np.zeros(3),'QUADRUPOLE':np.array([[0 ,0 ,0]]*3)}}", #,
            f"if do_mom:",
            ind+f"start_mom=time.time()",
            ind+f"mom= eval_moments(wfn)" ,
            *[ ind+x for x in store_timing(job, 'eval_mom', 'start_mom', add=True)],
            f"else: mom=None",
            f"{loc_store}.update(mom=mom)",
            # f"{mom_registry}.update({{",
            # ind+f"\'{pot_tag}\': mom",
            # "}})",
            f"return {loc_store},wfn"
    ]

    lines=[ 
        f"def get_moments({ext_pot_int}, {delta}, do_mom=True):",
        *[ ind+l for l in lines],
    ]
    return lines
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
        ind+f"fac=1",# if len(x)==4 else 2",
        ind+f"ar=-1*np.array(x).reshape(tuple([3]*int(np.emath.logn(3,x.size))))", # e-
        ind+f"return {tag_matrix.__name__}(ar) ",
        f"new_di= dict()",
        f"[ new_di.update(my_trans(wfn.variable(prop) ) ) for prop in {[ trans_prop(x) for x in properties]} ]",
        f"return new_di",
        # "return dict({ ",
        # *[ind+f"'{trans_prop(prop)}' :  my_trans(wfn.variable('{trans_prop(prop).upper()}')) ," for prop in properties], # the -1 is for the conversion from psi4 convention to usual convention (e.g. dipole is negative since it e-)
        # '})'
    ]

    lines=[ 
        f"def eval_moments({wfn_internal}):",
        *[ind+l for l in lines],
    ]
    return lines

def tensor_tag_from_indices(indices):
    """Return the tensor tag string for a given index combination (e.g. [0,1,2] -> 'xyz')."""
    mapper = ["x", "y", "z"]
    return ''.join(mapper[i] for i in indices)
def tag_matrix(M, field=None):
    def get_val(ind):
        return float(M[tuple(ind)])
    new_M={}
    M = np.asarray(M)
    base_tag=f"@{field}" if not field is None else ''

    dims=M.shape
    # get all indices of M
    all_indices = list(np.ndindex(M.shape))
    for ind in all_indices:
        tag=tensor_tag_from_indices(ind)
        new_M.update({tag+base_tag: get_val(ind)})

    return new_M


def pol_transform(pols):
    new_pols={}
    for t,pol in pols.items():
        for k,v in pol.items():
            new_tags=tag_matrix(v, t)
            if new_tags is None: continue
            for kk,vv in new_tags.items():
                if kk in new_pols.keys(): raise Exception(f"Already defined: \'{k}\'")
                else: new_pols.update({kk:vv})
    return new_pols

# Make a function that transform

def make_moments(source,the_var):
    return [
        f"# Get Multipoles from zero field calculation",
        f"{the_var}={{}}",
        f"for k,v in {source}.items():",
        f"    {the_var}.update({str(tag_matrix.__name__)}(v))",
    ]