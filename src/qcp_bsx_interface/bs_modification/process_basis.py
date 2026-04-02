from qcp_bsx_interface.util.import_helper import *

from .augmentation import (
    produce_augmentation, augment_basis, HIGHER_AUG_REGEX, AUGMENTATION_MAP
)
from .core_polarization import elements_require_core_polarization
from .util import recognized_basis, KNOWN_PARTICLES, formatted_basis_dict



@val_call
def bs_requires_modification(the_basis: recognized_basis) -> bool:
    """ Tell if basis name can be processed by psi4 or it needs to be explicitly be provided """
    found=[]
    for regex in [HIGHER_AUG_REGEX]+KNOWN_PARTICLES['conditional_core_polarization']:
        f=re.search(regex, the_basis.lower())
        if bool(f):
            found.append(f.group())
    return bool( len(found)>0 )

@val_call
def modify_bs(basis_name: recognized_basis, elements: List[str|int]) -> Dict:
    """ Expands the basis into a set of elements """

    # Return basis set of elements with basis set expanded (core polarized vs not core polarized )
    el_to_basis=elements_require_core_polarization(basis_name, elements)
    # Reverse that basis to list of elements
    basis_to_els=dict()
    for el,bs in el_to_basis.items():
        basis_to_els[bs]=basis_to_els.get(bs,[])+[el]

    bss={}
    els_there=set({})
    for bs, els in basis_to_els.items():
        if any( el in els_there for el in els ):
            raise ValueError(f"Element {el} appears in multiple basis sets, which is not supported (each element can only be augmented once, so if you want to augment an element with multiple basis sets, you need to call produce_augmentation multiple times).")
        els_there.add(frozenset(els))

        # If augmentation is required generate the bsx dictionary for this basis.
        if re.match(HIGHER_AUG_REGEX, bs.lower()):
            bss[bs]=produce_augmentation(bs,els)
        else:
            bss[bs]=els

    return bss


@val_call
def make_psi4_block(bases: Dict[str,list|Dict[str,Any]], basis_name='custom', is_df_basis:bool=False) -> List[str]:
    """
    The df_basis surpresses printing basis which does create a new basis which will not be used instead of changing the 
    df_basis!"""

    ind=4*' '
    

    def el_form(el):
        """ Format element (may be provided as int but needs to be str)"""
        if re.match(r'^\d+$', str(el)):
            el=nuclear_charge_to_element_symbol(int(el))
        if isinstance(el, int):
            el=nuclear_charge_to_element_symbol(el)
        return el
    def make_basis_str(basis_dict):
        new_bs_str_base=get_str(basis_dict,'psi4')
        #for atty,b in parts.items():
        return [ ind+ind+f"{x}" for x in new_bs_str_base.splitlines() if not x.strip()=='']
    def clarify_custom_bs(di):
        # make clear if basis sets are provided by prepending "custom_" tag
        bases_work=dict()
        for k,v in bases.items():
            if isinstance(v, dict):
                bases_work[f"custom_{k.replace('custom_','')}"]=v
            else:
                bases_work[k]=v
        return bases_work
    def get_dict_of_elements(basis_dict):
        el_to_basis={}
        for k,v in basis_dict.items():
            if isinstance(v, list):
                for el in v:
                    assert not el in el_to_basis, f"Element {el} appears in multiple basis sets, which is not supported (each element can only be augmented once, so if you want to augment an element with multiple basis sets, you need to call produce_augmentation multiple times)."
                    el_to_basis[el_form(el)]=k
            else:
                for el in v['elements'].keys():
                    assert not el in el_to_basis, f"Element {el} appears in multiple basis sets, which is not supported (each element can only be augmented once, so if you want to augment an element with multiple basis sets, you need to call produce_augmentation multiple times)."
                    el_to_basis[el_form(el)]=k
        return el_to_basis
    def custom_basis_set_declaration_lines(bs_di):
        custom_bs_lines=[]
        for bs,v in bs_di.items():
            if isinstance(v, dict):
                custom_bs_lines+=[f"[{bs}]"]
                custom_bs_lines+=make_basis_str(v)
        return custom_bs_lines
    
    # add custom prefix for custom basis sets
    bases_work=clarify_custom_bs(bases)
    el_to_basis=get_dict_of_elements(bases_work)
    # Check that all basis sets are covered by elements
    assert len(el_to_basis)>=len(bases), f"Number of elements covered by basis sets ({len(el_to_basis)}) is less than number of basis sets provided ({len(bases)}), which is not expected since each basis set should cover at least one element. Please check the input basis sets and make sure they are correctly formatted and that each basis set covers at least one element."

    # This will be the explicit declaration of basis sets
    custom_declaration=custom_basis_set_declaration_lines(bases_work)
    assignments=[
        (f"basis " if not is_df_basis else "")+f"{basis_name.replace('-','_')} {{", # identifier of the compounded basis sets
        # Assignment of basis sets to each element
        *[ ind + f"assign {nuclear_charge_to_element_symbol(el) if isinstance(el, int) else el} {bs}" for el,bs in el_to_basis.items()],
        # Explicit declaration of custom basis sets (if any)
        *[ ind+x for x in custom_declaration],
        f"}}",
    ]
    
    return assignments