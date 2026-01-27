from . import *

from qcp_objects.objects.properties import geometry

def psi4_molecule_string(atom_types, coordinates, units=None, total_charge=0, multiplicity=1, symm:bool=False):

    if isinstance(units, type(None)):
        units={'LENGTH':'BOHR'}
    content=[]

    # Comment (makes it easier to copy this block into xyz
    content+=[f"# Number of nuclei: {len(coordinates)}"]

    # Charge and multiplicity
    content+=[f"{total_charge} {multiplicity} # Charge and Multiplicity"]
    
    # Make geometry
    geometry=[]
    for at_ty, coord in zip(atom_types, coordinates):
        assert len(at_ty)<3 and at_ty!=''
        coord=[ f"{c:.8f}" for c in coord]
        coord=[ f"{c:>15}" for c in coord]
        geometry+=[ f"{at_ty.capitalize():<2} {' '.join(coord)}"]
    content+=geometry


    if not symm:
        content+=[ f"symmetry c1"]
    content+=[
        'nocom',
        'noreorient',
        f"units {''.join(list(units.values()))}"
    ]

    return '\n'.join(content)
class geometry_model(controller_model):
    """ Geometry class for psi4 jobs """
    ac_shift: float|None=None
    charge: int|None=None
    multiplicity: int|None=None
    xyz_file: str
    __example_values__=dict(
        xyz_file='molecule.xyz',
    )
    def check_geometry(self):
        if self.charge is None: raise Exception(f"Specify total charge (for safety)")
        if self.multiplicity is None: raise Exception(f"Specify multiplicity (for safety)")
    @val_call
    def molecule_lines(self, molecule_name:str, symm:bool=False) -> Tuple[List[str],str]:
        # self.units.LENGTH='ANGSTROM' # Import so that coordinates are returned in angstroms!
        geom=geometry(self.xyz_file)
        molecule_str='\n'.join([4*' '+x for x in psi4_molecule_string(
            atom_types=geom.atom_types, coordinates=geom.coordinates, total_charge=self.charge, multiplicity=self.multiplicity,
              units={'LENGTH':'ANGSTROM'}, symm=symm).split('\n') ],
        )

        assert molecule_name!='', "Provide a valid molecule name! (molecule name is empty)"
        mol_string=f"monomer_string"
        molecule_lines=[ f"#Define the molecule (definition in molecule string has the advantage of being able to create new molecules reliable)",
            f"{mol_string}=\"\"\"\n{molecule_str}\n\"\"\"", f"{molecule_name}=Molecule.from_string(monomer_string)"]
        return molecule_lines, mol_string
    @property
    def atom_types(self):
        geom=geometry(self.xyz_file)
        return geom.atom_types
    @property
    def coordinates(self):
        geom=geometry(self.xyz_file)
        return geom.coordinates
    @property
    def open_shell(self):
        """ Assumes that multiplicity 1 assumes open shell"""
        return (not self.multiplicity==1) 
    def __init__(self, **data):

        super().__init__(**data)
        try:
            geom=geometry(self.xyz_file)
            # self.atom_types=geom.atom_types
            # self.coordinates=geom.coordinates
        except Exception as e:
            raise ValueError(f"Invalid geometry file \'{self.xyz_file}\': {e}")
        self.check_geometry()

def wrap_geometry(v):
    if isinstance(v, geometry):
        ret = v
    elif isinstance(v, str):
        ret = geometry(geom=v)
    elif isinstance(v, dict):
        ret = geometry_model(**v)
    else:
        raise ValueError(f"Cannot wrap geometry from {v} of type {type(v)}")
    kwargs=ret.model_dump(include={'charge','multiplicity','ac_shift','xyz_file'})
    ret=geometry_model(
        **kwargs
    )
    ret.check_geometry()
    return ret
geometry_wrapper=Annotated[geometry_model, BeforeValidator(wrap_geometry)]
