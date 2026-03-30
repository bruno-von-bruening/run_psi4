from qcp_bsx_interface.util.import_helper import *
from .util import formatted_basis_dict, recognized_basis, KNOWN_PARTICLES

@val_call
def elements_require_core_polarization(
    basis_name: recognized_basis, elements: List[Union[str,int]]
) -> Dict[str,str]:
    

    # match the conditional group
    conditional_core_pol=None
    for regex in KNOWN_PARTICLES['conditional_core_polarization']:
        match=re.search(regex, basis_name.lower())  # re.search not re.match: pattern can appear anywhere in the string
        if bool(match):
            conditional_core_pol=match.group()

    base_name=basis_name.replace(conditional_core_pol,'') if conditional_core_pol else basis_name
    base_name=recognized_basis(base_name) # validate the base name (without the conditional core pol pattern) to make sure it's a recognized basis set, this also serves to validate the conditional core pol pattern since if it was not in the expected format, the base name would not be recognized after removing it
    base_core_polarized=re.match(r'cc-pv(wc|c)+[d,t,q,5,6]z', base_name.lower()) is not None
    
    elements=[ nuclear_charge_to_element_symbol(el) if isinstance(el, int) else el for el in elements ] 
    
    if conditional_core_pol is None:
        return { el:base_name for el in elements }
    else:
        assert not bool(base_core_polarized), f"Basis {basis_name} has conditional core polarization but base is already core polarized!"

        conditional_core_pol=conditional_core_pol.lstrip('-') # remove leading dash if present (since some patterns have it and some not, but we only care about the core polarization type and row number which are both after the dash, so we can remove it for simplicity)
        if conditional_core_pol.lower().startswith('wc'):
            type_of_core_pol='wc'
        elif conditional_core_pol.lower().startswith('c'):
            type_of_core_pol='c'
        else:
            raise Exception(f"Unexpected conditional core polarization pattern {conditional_core_pol} in basis name {basis_name}, expected format like wc3 or c3 where the number indicates the row for which to apply core polarization")

        # Get where to apply core polarization from the number in the pattern
        number_of_row=re.search(r'\dr', conditional_core_pol).group()
        try:
            number_of_row=int(number_of_row[:-1])
        except:
            raise Exception(f"Error parsing conditional core polarization pattern {conditional_core_pol} in basis name {basis_name}, expected format like wc3 or c3 where the number indicates the row for which to apply core polarization")
        cutoff_z={ # when each row ends (including index)
            1:0,
            2:2,
            3:10,
            4:18,
            5:36,
        }
        assert number_of_row in cutoff_z, f"Unexpected row number {number_of_row} parsed from conditional core polarization pattern {conditional_core_pol} in basis name {basis_name}, expected an integer between 1 and 5 indicating the row for which to apply core polarization (e.g., wc3 or c3 for third row elements)"
        cutoff_z=cutoff_z[number_of_row]
        Zs=[ element_symbol_to_nuclear_charge(el) for el in elements ]
        basis_di={}
        for Z,el in zip(Zs,elements):
            if Z <= cutoff_z:
                basis_di[el]=base_name
            else:
                the_name=re.sub(r'cc-pv', f'cc-p{type_of_core_pol}v', base_name, flags=re.IGNORECASE)
                the_name=recognized_basis(the_name) # validate the modified basis name to make sure it's recognized, this also serves to validate the conditional core pol pattern since if it was not in the expected format, the modified name would not be recognized after adding it in)
                basis_di[el]=the_name
        return basis_di