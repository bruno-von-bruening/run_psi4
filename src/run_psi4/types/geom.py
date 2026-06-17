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
    charge: int
    multiplicity: int
    xyz_file: str
    _geometry: geometry|None=None
    __example_values__=dict(
    )
    @classmethod
    def is_dummy(cls, input:str):
        return input.lower().startswith('dummy')
    def check_geometry(self):
        if self.charge is None: raise Exception(f"Specify total charge (for safety)")
        if self.multiplicity is None: raise Exception(f"Specify multiplicity (for safety)")
    @val_call
    def molecule_lines(self, molecule_name:str, symm:bool=False) -> Tuple[List[str],str]:
        # self.units.LENGTH='ANGSTROM' # Import so that coordinates are returned in angstroms!
        geom=self._geometry
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
        return self._geometry.atom_types
    @property
    def coordinates(self):
        return self._geometry.coordinates
    @property
    def open_shell(self):
        """ Assumes that multiplicity 1 assumes open shell"""
        return (not self.multiplicity==1) 
    def __init__(self, **data):

        super().__init__(**data)
        if self.is_dummy(self.xyz_file):
            try:
                matcher= r"dummy_[a-zA-Z1-9\-]+"
                assert bool(re.fullmatch(matcher, self.xyz_file)), f"Provided input startswith dummy but does not match expected format for dummy values: {matcher}"

                dummy_string, kind=self.xyz_file.split('_')
                avail_dummy_kinds=['h2','be','he','ne','ar','h2o']
                if kind.lower() not in avail_dummy_kinds:
                    raise ValueError(f"Unrecognized dummy kind '{kind}'. Available dummy kinds are: {avail_dummy_kinds}")

                if kind.lower()=='h2':
                    geom=geometry(atom_types=['H','H'], coordinates=np.array([ 0]*5+[.74]).reshape(-1,3), length_units='ANGSTROM')
                elif kind.lower() in ['be','ne','ar','he']:
                    geom=geometry(atom_types=[kind.capitalize()], coordinates=np.array([[0,0,0]]), length_units='ANGSTROM')
                elif kind.lower()=='h2o':
                    xyz="""O 	0.0000 	0.0000 	0.1173
                    H 	0.0000 	0.7572 	-0.4692
                    H 	0.0000 	-0.7572 	-0.4692"""
                    atom_types,coordinates=zip( *[ (x.split()[0], [ float(x) for x in x.split()[1:]]) for x in xyz.strip().split('\n')])
                    geom=geometry(atom_types=atom_types, coordinates=np.array(coordinates, dtype=float), length_units='ANGSTROM')
                else:
                    raise Exception(f"Unrecognized dummy tag '{kind}'")
            except Exception as ex: raise Exception(f"Could not interpret {self.xyz_file}: {ex}") from ex
        else:
            try:
                geom=geometry(self.xyz_file)
            except Exception as e:
                raise ValueError(f"Invalid geometry file \'{self.xyz_file}\': {e}")
            # self.atom_types=geom.atom_types
            # self.coordinates=geom.coordinates
        self._geometry=geom
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
