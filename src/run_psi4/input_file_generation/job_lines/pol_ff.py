from . import *

from run_psi4.types.job import psi4_job



@val_call
def make_external_field(tensor_key:str, delta:float) -> Tuple[List[str], str,str]:
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
    lines+=[f"psi4.set_options({{",
     f"\"scf_type\": \"df\",",
    f"\"e_convergence\": 1e-12,",      # energy convergence
    f"\"d_convergence\": 1e-12,",       # density matrix convergence
    f"}})"]
    lines+=[f"{delta_key}={delta}", f"{d_key}={d}"]
    if len(tensor_key)==1:
        vec=get_dir(tensor_key)
        lines+=[
            f"vec=np.array([ {vec[0]}, {vec[1]}, {vec[2]}])",
            f"vec_inf={d_key}*vec",
            f"q={d_key}**2/2 * {delta_key} ",
        ]
        vecs=[
            ( -1,"-vec_inf"),
            ( 1,"vec_inf"),
        ]

    elif len(tensor_key)==2:
        the_t=tensor_key[0]
        vec=get_dir(tensor_key[0])
        sep=1.e-2
        lines+=[
            f"vec=np.array([ {vec[0]}, {vec[1]}, {vec[2]}])",
            f"print_out(f\"COM is: {{the_molecule.center_of_mass()}}\")",
            f"sep={sep}",
            f"x=np.array([ 1. ,0. ,0. ])",
            f"y=np.array([ 0. ,1. ,0, ])",
            f"z=np.array([ 0. ,0. ,1. ])",
            # f"vec_inf={d_key}*vec",
            # f"vec_mux=sep*x",
            # f"vec_muz=sep*z",
            # f"vec_muy=sep*y",
            f"q=-({d_key}**3)/(8) * {delta_key} *2   ",
            f"print_out( f\"q is {{q}}\")",
            f"coords=the_molecule.full_geometry().to_array()",
            f"com=np.sum([ the_molecule.Z(i)/np.sum([the_molecule.Z(i) for i in range(len(coords))]) * c for i,c in enumerate(coords) ] ,axis=0)",
            #f"q=({d_key}**4)/(8*sep) * {delta_key}  ",
            #f"center=np.array([ 1, 0 ,0 ])",
            f"center=np.array([ 0, 0 ,0 ])",
            f"center=the_molecule.center_of_mass()",
            f"center=np.array([ center[0],center[1],center[2]])",
            f"center=com",
            f"print_out(f\"COM is {{com}} \")",
            #f"q=1/(6*{d**2) * {delta_key}  "
        ]
        vecs_bare=[]
        if True:
            for t in []:#['x']:#,'y','z']:
                vecs_bare+=[
                    # (-.5, f"center + {t}*d     - {t}*sep/2"),
                    # (+.5, f"center + {t}*d     + {t}*sep/2"),
                    # (+.5, f"center - {t}*d    + {t}*sep/2"),
                    # (-.5, f"center - {t}*d    - {t}*sep/2"),
                    (1., f"center - {t}*d    "),
                    (1., f"center + {t}*d    "),
                    # (-.5, "vec_inf"),
                    # (+.5, "-vec_inf"),
                    # (+.5* 1, "vec_inf+vec_mu"),
                    # (-.5* 1, "-vec_inf-vec_mu"),
                    # (1, "d*x"),
                    # (1, "-d*x"),
                    # (-.5, "-d*y"),
                    # (-.5, "+d*y"),
                    # (-.5, "-d*z"),
                    # (-.5, "+d*z"),
                ]
        if True:
            for t in ['y','z','x']:
                if t==the_t:
                    continue
                fac=-1
                vecs_bare+=[
                    # (-.5*fac, f"center + {t}*d     - {t}*sep/2"),
                    # (+.5*fac, f"center + {t}*d     + {t}*sep/2"),
                    # (+.5*fac, f"center - {t}*d    + {t}*sep/2"),
                    # (-.5*fac, f"center - {t}*d    - {t}*sep/2"),
                    (1*fac, f"center - {t}*d  "),
                    (1*fac, f"center + {t}*d  "),
                ]
        vecs=[]
        #vecs=vecs_bare
        #tensor_remainder
        #for x in ['x','y','z']:
        #    if not x== tensor_key
        #for q,v in vecs_bare:
        #    vecs+=[ (q,v+x) for x in ['+z','-z','+y','-y']]
        vecs=vecs_bare
        print(vecs_bare)
    
    pot_key_pl=f"external_potential_{tensor_key}_pl"
    pot_key_mn=f"external_potential_{tensor_key}_mn"

    ind=' '*4
    lines+=[
        f"pot_base=[",
    ]+[
        ind+f"({x}, {tag})," for x,tag in vecs 
    ]+[
        ']',
        f"{pot_key_pl}=[ (x*(q),v  ) for x,v in pot_base ] ",
        f"{pot_key_mn}=[ (x*(-q),v ) for x,v in pot_base ] ",
    ]
    var_mem.update({
        'external_pot_mn':pot_key_mn,
        'external_pot_pl':pot_key_pl,
        'delta':delta_key,
    })

    return lines

def get_dipole_lines(job):
    lines=[]
    if True:
        my_kwargs=dict(
            molecule=job.var_names.molecule_name,
            return_wfn=True,
        )

        op='gradient'
        lines+=[ f"_, wfn = {op}({job.method.func_var}, "+' , '.join( f"{k}={v}" for k,v in my_kwargs.items() )+")"
        ]


        tag='QUADRUPOLE'
        #tag='DIPOLE'
        dp_var=f"ind_dp_zero"
        var_mem.update( {dp_var:dp_var})
        lines+=[
            f"oeprop = psi4.core.OEProp(wfn)",
            f"oeprop.add('{tag}')",
            f"oeprop.compute()",
            f"{dp_var}=psi4.core.variable('{tag}')",
            f"psi4.core.clean()",
        ]
    for dir in ['mn','pl']:
        my_kwargs=dict(
            molecule=job.var_names.molecule_name,
            external_potentials=get_var(f"external_pot_{dir}"),
            return_wfn=True,
        )

        op='gradient'
        lines+=[ f"_, wfn = {op}({job.method.func_var}, "+' , '.join( f"{k}={v}" for k,v in my_kwargs.items() )+")"
        ]


        tag='QUADRUPOLE'
        #tag='DIPOLE'
        dp_var=f"ind_dp_{dir}"
        var_mem.update( {dp_var:dp_var})
        lines+=[
            f"oeprop = psi4.core.OEProp(wfn)",
            f"oeprop.add('{tag}')",
            f"oeprop.compute()",
            f"{dp_var}=psi4.core.variable('{tag}')",
            f"psi4.core.clean()",
     ]

    lines+= [ 
        f"print_out( 'plus\\n' )", 
        f"print_out( str({get_var('ind_dp_pl')})+'\\n' )", 
        f"print_out( 'minus\\n' )", 
        f"print_out( str({get_var('ind_dp_mn')}) )",
        f"print_out( 'zero\\n' )", 
        f"print_out( str({get_var('ind_dp_zero')}) )",
        f"pol= ( {get_var('ind_dp_pl')}-{get_var('ind_dp_zero')} )/ ({get_var('delta')}) ",
    ]
    

    return lines


def generate_pol_ff_job_lines(job:psi4_job)->list[str]:
    """Generate job lines to compute polarisability tensors and store them in a yaml file.

    Args:
        job (psi4_job): The psi4 job object.
    """

    lines=[]

    # for plus minus

    #for tag in ['x','y','z', 'xx', 'yy', 'zz']:
    for tag in ['xx','yy','zz']:
        
        # Define reference points and memorise the results

        pot_lines =make_external_field(tag,1.0e-10) # makes external corresponding to perturbation

        dipole_lines=get_dipole_lines(job)

        lines+=pot_lines+dipole_lines
        lines+=[ f"print_out( f\"\\n{tag}: pol is \\n"+
                '\\n'.join([
                    f"{{pol[{i}][0]:.6f}} {{pol[{i}][1]:.6f}} {{pol[{i}][2]:.6f}}" for i in range(3) 
                ])
                +"\\n\")"
        ]



    # Perform single point within field

    # Get the dipole difference

    return lines

    # Generate finite field tensor (should auto generate this one)


    # Mkae hamiltonian
    # set_options( {'perturb_dipole': make_c2v_vec(delta),'reference':'rhf' })
    # set {
    # OPDM_RELAX True
    # }

    # Single point

    # Get dipole

    # Make finite field