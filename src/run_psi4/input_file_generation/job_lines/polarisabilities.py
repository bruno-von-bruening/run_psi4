from . import *

from run_psi4.types.job import psi4_job

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
    f"alpha = wfn.variable('CCSD DIPOLE POLARIZABILITY')",
    f"quit()",
    ]

    # Get the polarisabilities
    def recover_polarisabilities(origin, type_of_origin=None, units_of_origin=None):
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
        expected_subindices={{
            'dipole':make_product(['x','y','z']),
            'quadrupole':make_product(['xx','yy','zz','xz','xy','yz']),
        }}

        polarisability_tensors={{}}
        zero_thres=1.e-12
        for rank, rank_key in key_per_rank.items():
        \ttensor={{}}
        \tavailable_keys=[x for x in the_variables.keys() if x.startswith(rank_key)]
        \tcovered_keys=dict([(x,False) for x in available_keys ])
        \tfor index in expected_subindices[rank]:
            \t\tfull_index=' '.join([rank_key,index.upper()])
            \t\tassert full_index in available_keys, [index, available_keys, rank_key, the_variables.keys()]
            \t\t# Mark the index as covered
            \t\tassert not covered_keys[full_index], f"Double read in of {{full_index}}?"
            \t\tcovered_keys[full_index]=True
            \t\t# Safe the index and value into the tensor
            \t\tvalue=the_variables[full_index]
            \t\tvalue=0 if abs(value)<zero_thres else value
            \t\ttensor.update({{index:value}})
        \tpolarisability_tensors.update({{rank:tensor}})
        \tassert all(covered_keys.values())

        # Drop to file
        polarisability_tensors.update({{ 'origin': {{'coordinate':{COM}, 'type_of_origin':\'{type_of_origin}\', 'units_of_origin':\'{units_of_origin}\' }}}})
        polarizability_file='polarisability_tensors.yaml'
        with open(polarizability_file,'w') as wr:
          \tyaml.dump(polarisability_tensors, wr)
        {job.var_names.storage}['files'].update({{'polarizability_tensors':polarizability_file}})
        """
        ind=4*' '
        recover_polarisabilities=[x.strip(' ') for x in recover_polarisabilities.split('\n')]
        return recover_polarisabilities
        # Sorting according to stones convention (See Cartesian–Spherical Conversion Tables@Stone's Theory of Intermolecular Interactions book)

    define_polarisability_origin=[x.strip() for x in f"""
    # Define center of mass (in angstrom!) and assert that it has been read in correctly!
    core.set_global_option("PROPERTIES_ORIGIN", {COM})
    assert {COM}==get_option('None','PROPERTIES_ORIGIN'), 'Check Center of mass'
    """.split('\n')]
    # Pierce file together
    lines=(        []
        + define_polarisability_origin
        + property_calculation
        + recover_polarisabilities(COM, type_of_origin='COM',units_of_origin=units_of_origin)
    )
    return lines