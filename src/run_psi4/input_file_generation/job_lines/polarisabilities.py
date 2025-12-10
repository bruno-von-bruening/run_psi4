from . import *

from run_psi4.types.job import psi4_job, store


def get_polarizability(job:psi4_job):
    """ Given a geometry calculate the Polarisabilities available in Psi4
    [ ] Provide starting guess (idk if that will be that much faster though)
    """

    # Get method information
    func=job.method.method
    basis_set=job.basis_set
    # Get geometry information
    molecule_name=job.var_names.molecule_name
    atom_types=job.geom.atom_types
    coordinates=job.geom.coordinates

    def calculate_center_of_mass(atom_types, coordinates):
        """ 
        Calculate center of mass as
        x_com = sum ( w_i * x_i )    with x_i beeing coordinates and w_i weights: w_i = Z_i/( sum Z_j )
        """
        atom_numbers=[]
        for at_ty in atom_types:
            if isinstance(at_ty, int):
                atom_numbers.append(at_ty)
            elif isinstance(at_ty, str):
                assert len(at_ty)<3 and len(at_ty)>0
                Z=element_symbol_to_nuclear_charge(at_ty)
                atom_numbers.append(Z)
        atom_numbers=np.array(atom_numbers)
        com=np.sum(atom_numbers[:,None]*coordinates, axis=0)/np.sum(atom_numbers)
        return com


    COM, units_of_origin=[ [ float(x) for x in calculate_center_of_mass(atom_types, coordinates)], 'ANGSTROM' ] # recast that for easy printout


    property_calculation=[
    "# perform polarizability calculation",
    "props = ['DIPOLE_POLARIZABILITIES','QUADRUPOLE_POLARIZABILITIES']",
    #"set_options({",
    #"   \'basis\' : basis_set,",
    #"   \'reference\' : 'rks',",
    #"})",
    #"print(method, monomer, props)",
    f"E, wfn = properties(method, molecule={job.var_names.molecule_name}, properties= props, return_wfn=True)",
    # f"alpha = wfn.variable('CCSD DIPOLE POLARIZABILITY')",
    # f"quit()",
    ]

    # Get the polarisabilities
    def recover_polarisabilities(origin, the_var:str, type_of_origin=None, units_of_origin=None):
        recover_polarisabilities=f"""
        #### Recover the polarisabilities saved in the global variables
        # Manipulate the 'variables' variable where polarisabilities are stored
        the_variables=variables()
        key_per_rank={{
            'dipole':'DIPOLE POLARIZABILITY',
            'quadrupole':'QUADRUPOLE POLARIZABILITY'
        }}
        if any( method.lower().startswith(x) for x in [ 'cc','mp'] ): key_per_rank=dict(
        (k,method.upper()+' '+v) for k,v in key_per_rank.items()
        )
        # Construct all expected subindices
        def make_product(the_list):
        \tnew_list=[]
        \tfor x in the_list:
            \t\tfor y in the_list:
            \t\t\tnew_list.append(x+y)
        \treturn new_list
        dip_entries=['x','y','z']
        quad_entries=['xx','yy','zz','xy','xz','yz']
        expected_subindices=[]
        for x in [dip_entries, quad_entries]:
            \texpected_subindices+=[
            \t\t  (a,b) for a in x for b in x
            \t]
        covered_keys=dict( (k,False) for k in expected_subindices )
        {the_var}=dict()
        for a,b in expected_subindices:
            \t#full_index=' '.join([rank_key,index.upper()])
            \t#assert full_index in available_keys, [index, available_keys, rank_key, the_variables.keys()]
            \tfull_index=' '.join([key_per_rank['dipole' if len(a)==1 else 'quadrupole'].upper(), a.upper()+b.upper()])
            \t# Mark the index as covered
            \tassert not covered_keys[(a,b)], f"Double read in of {{full_index}}?"
            \tcovered_keys[(a,b)]=True
            \t# Safe the index and value into the tensor
            \tvalue=the_variables[full_index]
            \t#value=0 if abs(value)<zero_thres else value
            \t{the_var}.update({{f"{{a}}@{{b}}":value}})
        assert all(covered_keys.values())"""
        return [x.lstrip(" ") for x in recover_polarisabilities.split('\n')]

    #for rank, rank_key in key_per_rank.items():
    #\ttensor={{}}
    #\tavailable_keys=[x for x in the_variables.keys() if x.startswith(rank_key)]
    #\tcovered_keys=dict([(x,False) for x in available_keys ])
    #\tfor index in expected_subindices[rank]:
    #    \t\tfull_index=' '.join([rank_key,index.upper()])
    #    \t\tassert full_index in available_keys, [index, available_keys, rank_key, the_variables.keys()]
    #    \t\t# Mark the index as covered
    #    \t\tassert not covered_keys[full_index], f"Double read in of {{full_index}}?"
    #    \t\tcovered_keys[full_index]=True
    #    \t\t# Safe the index and value into the tensor
    #    \t\tvalue=the_variables[full_index]
    #    \t\tvalue=0 if abs(value)<zero_thres else value
    #    \t\ttensor.update({{index:value}})
    #\t{the_var}.update({{rank:tensor}})
    #\tassert all(covered_keys.values())

#         # Drop to file
#         polarisability_tensors.update({{ 'origin': {{'coordinate':{COM}, 'type_of_origin':\'{type_of_origin}\', 'units_of_origin':\'{units_of_origin}\' }}}})
#         polarizability_file='polarisability_tensors.yaml'
#         with open(polarizability_file,'w') as wr:
#           \tyaml.dump(polarisability_tensors, wr)
#         {job.var_names.storage}['files'].update({{'polarizability_tensors':polarizability_file}})
#         """
#         ind=4*' '
#         recover_polarisabilities=[x.strip(' ') for x in recover_polarisabilities.split('\n')]
#         return recover_polarisabilities
#         # Sorting according to stones convention (See Cartesian–Spherical Conversion Tables@Stone's Theory of Intermolecular Interactions book)

    from .pol_ff_functions import center_of_charge_function
    # Get center of charge to define multipoles around
    func_name, new_lines = center_of_charge_function()
    lines+=new_lines
    center_of_charge_var='center_of_charge'
    var_mem.update({'center_of_charge':center_of_charge_var})
    lines+=[
        f"{get_var('center_of_charge')} = {func_name}({job.var_names.molecule_name})",
    ]
    define_polarisability_origin=new_lines +[x.strip() for x in f"""
    # Define center of mass (in angstrom!) and assert that it has been read in correctly!

    core.set_global_option("PROPERTIES_ORIGIN", {get_var('center_of_charge')})
    assert {COM}==get_option('None','PROPERTIES_ORIGIN'), 'Check Center of mass'
    """.split('\n')]
    # Pierce file together
    lines=(        []
        + define_polarisability_origin
        + property_calculation
        + recover_polarisabilities(COM, 'pol_linresp',type_of_origin='COM',units_of_origin=units_of_origin)
        + store(job, 'polarisabilities_linear_response', 'pol_linresp')
    )
    return lines