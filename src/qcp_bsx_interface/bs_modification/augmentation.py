from qcp_bsx_interface.util.import_helper import *

from .util import formatted_basis_dict, HIGHER_AUG_REGEX, recognized_basis, AUGMENTATION_MAP

def augment_basis(basis_name: str, el: str | int):
    # First character indicates augmentation level (has to be [d,t,q]-aug-) as checked by requires_augmentation
    level = basis_name[0]
    aug_level = AUGMENTATION_MAP[level]
    base = basis_name.replace(f'{level}-aug-', 'aug-')

    # For H the aug-cc-pc basis is not available (since there is no core) need to fix that here
    Z = element_symbol_to_nuclear_charge(el) if not isinstance(el, int) else el
    if Z == 1 and bool(re.match(r'.*cc-p[w]?cv.*', basis_name.lower())):
        raise Exception("Augmented basis for H is not available (Should provide the regular version beforehand!)")

    old_bs_dict = get_basis(base, elements=[el])
    return geom_aug(old_bs_dict, aug_level, use_copy=True)

@val_call
def produce_augmentation(basis_name: List[str] | str, elements: List[str | int]) -> formatted_basis_dict:

    # Get an list of unique pairs between bases names and elements
    if isinstance(basis_name, str):
        basis_name=[basis_name for _ in elements]
    assert len(basis_name)==len(elements), f"Length of basis_name list ({len(basis_name)}) must match length of elements list ({len(elements)})"
    bs_el_pairs=list(set([ (bs,el) for bs,el in zip(basis_name, elements)  ])  )

    # 
    bss=dict()
    for bs,el in bs_el_pairs:
        # At this stage an element reoccurring means it has two different basis sets assigned
        if el in bss: 
            raise ValueError(f"Element {el} appears multiple times in the input with different basis sets, which is not supported (each element can only be augmented once, so if you want to augment an element with multiple basis sets, you need to call produce_augmentation multiple times).")

        # Provide standard or augmented basis set
        highaug=re.match(HIGHER_AUG_REGEX, bs.lower())
        if bool(highaug):
            bss[el] = augment_basis(bs, el)
        else:
            bss[el] = get_basis(bs, elements=[el])

    # Now combine them
    assert len(bss)>0, f"Did not produce any basis set augmentation for basis {basis_name} and elements {elements}"

    for i,(k,v) in enumerate(bss.items()):
        if i==0:
            bs_dict=v # take meta data from first entry
        else:
            for atty,el_data in v['elements'].items():
                assert not atty in bs_dict['elements'], f"Unexpectedly found the same element {atty} in multiple basis sets (this should not happen since we are only augmenting one element at a time, but it seems to have happened here, so I will raise an error to avoid silently producing wrong results)"
                bs_dict['elements'][atty]=el_data

    return bs_dict