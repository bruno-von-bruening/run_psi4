from . import *

# def make_multipolar_field():


@val_call
def make_external_field(tensor_key:str) -> List[str]:
    """ Most basic unfiltered or oriented tensors for given rank """
    lines=[]

    _max_rank=2 # only quadrupoles supported for now
    _min_rank=0 # ALlowd monopoles (may need that for some tests of linear response)

    if len(tensor_key)==0:
        raise ValueError("Tensor key \'{tensor_key}\' must have at least one element indicating the rank of the tensor.")
    elif len(tensor_key)<_min_rank:
        raise ValueError(f"Tensor key \'{tensor_key}\' length is below minimum supported rank of {_min_rank}.")
    elif len(tensor_key)>_max_rank:
        raise ValueError(f"Tensor key \'{tensor_key}\' length exceeds maximum supported rank of {_max_rank}.")
    else:
        pass

    def get_dir(tensor_key:str):
        assert len(tensor_key)==1, f"Expected input as single character"
        mapper=dict( (v,i) for i,v in enumerate(['x','y','z']) )
        if not tensor_key in mapper:
            raise ValueError(f"Tensor key \'{tensor_key}\' is not a valid key.")
        pos=mapper[tensor_key]
        vec=np.zeros(3) ; vec[pos]=1
        return vec

    cutoff=10 # bohr
    d=1000# 1% error margin

    vec_mu=None
    delta_key='delta'; d_key='d'
    lines+=[
        f"coords=the_molecule.full_geometry().to_array()",
        f"center=np.sum([ the_molecule.Z(i)/np.sum([the_molecule.Z(i) for i in range(len(coords))]) * c for i,c in enumerate(coords) ] ,axis=0)",
        f"x=np.array([ 1. ,0. ,0. ])",
        f"y=np.array([ 0. ,1. ,0, ])",
        f"z=np.array([ 0. ,0. ,1. ])",
    ]
    lines+=[ f"{d_key}={d}"]
    if len(tensor_key)==1:
        vec=get_dir(tensor_key)
        lines+=[
            f"vec=np.array([ {vec[0]}, {vec[1]}, {vec[2]}])",
            f"vec_inf={d_key}*vec",
            f"q=-{d_key}**2/2 ",#* {delta_key} ",
        ]
        vecs_bare=[
            ( -1, f"center - d*{tensor_key}"),
            ( 1,f"center + d*{tensor_key}"),
        ]

    elif len(tensor_key)==2:

        vec=get_dir(tensor_key[0])

        vecs_bare=[]
        a,b=tensor_key[0],tensor_key[1]

        if a==b:
            fac=3
        else:
            fac=6
        lines+=[
            f"q=-({d_key}**3)/({fac})     ",
        ]

        if a==b:
            vecs_bare+=[
                (1., f"center - {a}*d    "),
                (1., f"center + {a}*d    "),
            ]
        else:
            fac=f"pref"
            lines+=[f"{fac}=d/np.sqrt(2)"]
            v_help=[                    # The vectors (1,1)^T and (1,-1)^T
                ["+", "+", +1],
                ["-", "+", -1 ]
            ]
            def flip(x):
                if x=='+': return '-'
                elif x=='-': return '+'
                else: raise Exception(x)
            for fi,se, sign in v_help:
                vecs_bare+=[
                    (1.*sign, f"center {fi} {fac}*{a} {se} {fac}*{b}"),
                    (1.*sign, f"center {flip(fi)} {fac}*{a} {flip(se)} {fac}*{b}")
                ]
    else: raise Exception()


    lines+=[
        f"pot_base=dict(",
        ind+f"charges = ["
    ] + [
        ind+ind+f"({x}*q, {tag})," for x,tag in vecs_bare  
    ]+[
        ind+f"],",
        ind+f"center = center,",
        ind+f"separation = d,",
        f")",
    ]
    lines+=[f"{get_var('pot_registry')}.update({{ '{tensor_key}':pot_base }})"]
    
    return lines

def get_moment_function(job):
    """ For each perturbation tag:
        1. get the right field
        2. compute the wfn
        3. evaluate moments
        4. store moments in registry
    """
    # first get the wfn
    lines=[]

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
    )
    action='gradient' # only needed if not scf or KS, otherwise gradient updates density matrix
    lines+=[ psi4_action(job, action, **kwargs)]

    # 3. evaluate moments
    lines+=[ 
            '',
            f"mom=eval_moments(wfn)" ,
            # f"{mom_registry}.update({{",
            # ind+f"\'{pot_tag}\': mom",
            # "}})",
            f"return mom"
    ]

    lines=[ 
        f"def get_moments({ext_pot_int}, {delta}):",
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
    properties=['DIPOLE','QUADRUPOLE']
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
        return tag
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
        f"return dict( ",
        *[ind+f"{prop} = psi4.core.variable('{trans_prop(prop).upper()}')," for prop in properties],
        ')'
    ]

    lines=[ 
        f"def eval_moments({wfn_internal}):",
        *[ind+l for l in lines],
    ]
    return lines

def pol_transform(pols):
    def tag_matrix(M, field):
        def get_val(ind):
            return float(M[tuple(ind)])
        new_M={}

        M = np.asarray(M)
        mapper=["x", "y", "z"]
        ndim=len(M.shape)
        if ndim==1:
            assert len(M)==3
            n=len(M)
            m=1
        elif ndim==2:
            # M=M-np.eye(M.shape[0])*np.trace(M)/3
            assert M.shape==(3,3)
            n, m = M.shape
        else: raise Exception(f"Unhandled dimension number: \"{ndim}\": {M}")

        base_tag=f"@{field}"
        for i in range(n):
            if ndim>1:
                for j in range( max(m,i+1)):
                    new_M.update({ mapper[i]+mapper[j]+base_tag:get_val([i,j])})
            else:
                new_M.update({ mapper[i]+base_tag: get_val([i])})
        return new_M


    new_pols={}
    for t,pol in pols.items():
        for k,v in pol.items():
            new_tags=tag_matrix(v, t)
            for kk,vv in new_tags.items():
                if kk in new_pols.keys(): raise Exception(f"Already defined: \'{k}\'")
                else: new_pols.update({kk:vv})
    return new_pols