from . import *
from qcp_bsx_interface.util.import_helper import *

from qcp_bsx_interface.augmentation.produce_augmentation import produce_augmentation, requires_augmentation



def count_basis_functions_for_element(element, basis_name):
    """
    Count number of basis functions for a single atom
    """
    if requires_augmentation(basis_name):
        basis_dict = produce_augmentation(basis_name, elements=[element])
    else:
        basis_dict = bsx.get_basis(basis_name, elements=[element])

    element_data = basis_dict['elements'][str(bsx.lut.element_Z_from_sym(element))]
    shells = element_data['electron_shells']

    print(bsx.write_formatted_basis_str(basis_dict, 'psi4'))

    nbf = 0
    for shell in shells:
        print(shell)
        if len(shell['exponents']) > 1:
            fac=len(shell['coefficients'])
        else:
            fac=1
        for l in shell['angular_momentum']:
            print(l)
            nbf += fac*(2 * l + 1)

    return nbf

def parse_formula(formula):
    """
    Parse a chemical formula like 'H2O' into a dict {'H': 2, 'O': 1}
    """
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    composition = dict()
    for element, count in tokens:
        composition[element] = composition.get(element, 0) + ( int(count) if count else 1 ) 
    return composition

def element_count(elements:dict[str,int]|List[str]|str) -> dict[str,int]:
    if isinstance(elements, str):
        elements=parse_formula(elements)
    if isinstance(elements, list): # count occurrences
        elem_count=dict()
        for el in elements:
            elem_count[el]=elem_count.get(el,0)+1
        elements=elem_count
    return elements

# Calculate the number of basis functions for a given molecule
@val_call
def calculate_nbasis(
        basis_set:str|dict, elements:dict[str,int]|List[str]|str
)->int:
    elements=element_count(elements)

    el_cost=dict()
    for el, num in elements.items():
        el_cost[el]=count_basis_functions_for_element(el, basis_set)

    return sum(el_cost.get(el, 0) * cnt for el, cnt in elements.items())
