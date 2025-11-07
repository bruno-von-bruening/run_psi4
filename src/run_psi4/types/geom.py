from . import *

from qcp_objects.objects.properties import geometry

def psi4_molecule_string(atom_types, coordinates, units=None, total_charge=0, multiplicity=1):

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


    content+=[
        'Symmetry C1',
        'nocom',
        'noreorient',
        f"units {''.join(list(units.values()))}"
    ]

    return '\n'.join(content)
class geometry(geometry):
    ac_shift: float|None=None
    def check_geometry(self):
        if self.charge is None: raise Exception(f"Specify total charge (for safety)")
        if self.multiplicity is None: raise Exception(f"Specify multiplicity (for safety)")
    @val_call
    def molecule_lines(self, molecule_name:str):
        self.units.LENGTH='ANGSTROM' # Import so that coordinates are returned in angstroms!
        molecule_str='\n'.join([4*' '+x for x in psi4_molecule_string(
            atom_types=self.atom_types, coordinates=self.coordinates, total_charge=self.charge, multiplicity=self.multiplicity, units={'LENGTH':'ANGSTROM'}).split('\n') ])
        assert molecule_name!='', "Provide a valid molecule name! (molecule name is empty)"
        molecule_lines=[ f"#Define the molecule (definition in molecule string has the advantage of being able to create new molecules reliable)",
            f"monomer_string=\"\"\"\n{molecule_str}\n\"\"\"", f"{molecule_name}=Molecule.from_string(monomer_string)"]
        return molecule_lines
    @property
    def open_shell(self):
        """ Assumes that multiplicity 1 assumes open shell"""
        return (not self.multiplicity==1) 

def wrap_geometry(v):
    if isinstance(v, geometry):
        ret = v
    elif isinstance(v, str):
        ret = geometry(geom=v)
    elif isinstance(v, dict):
        ret = geometry(**v)
    else:
        raise ValueError(f"Cannot wrap geometry from {v} of type {type(v)}")
    ret.check_geometry()
    return ret
geometry_wrapper=Annotated[geometry, BeforeValidator(wrap_geometry)]