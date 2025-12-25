from . import *

def center_of_charge_function():
    """ Return the function name and lines defining it
    The function computes the center of charge of a molecule
    """
    lines=[]
    func_name='get_center_of_charge'
    lines+=[
        f"def {func_name}(mol):",
        ind+f"coords=mol.full_geometry().to_array()",
        ind+f"center=np.sum([ mol.Z(i)/np.sum([mol.Z(i) for i in range(len(coords))]) * c for i,c in enumerate(coords) ] ,axis=0)",
        ind+f"return center",
    ]
    return func_name, lines
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