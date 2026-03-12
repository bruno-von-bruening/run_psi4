
from basis_set_exchange.manip import geometric_augmentation as geom_aug
from basis_set_exchange  import write_formatted_basis_str as get_str, get_basis

from typing import List
import re
augmentation_map = {
    'd': 1,
    't': 2,
    'q': 3
}

def requires_augmentation(the_basis: str) -> bool:
    return bool(re.match(r'^[d,t,q]-aug-', the_basis.lower()))


def modify_basis(basis_dict):
    import copy
    new_bas=copy.deepcopy(basis_dict)
    elements=new_bas['elements']
    for atty,el_data in elements.items():
        vals=el_data['electron_shells']
        new_shells=[]

        # need to remeber last
        last_exp=None
        for i,shell in enumerate(vals):


            if 'exponents' in shell:
                # If the new exponent repeats the last contracted then delete this contraction!
                shell['exponents']=[ f"{round(float(x),5):.7f}" for x in shell['exponents'] ]
                    
            if 'coefficients' in shell:
                if len(shell['coefficients'])>1:
                    #assert len(shell['coefficients'])==2, f"Unexpected number of contractions in shell (expected 2 one for flexible diffuse, one for core): {len(shell['coefficients'])}: {shell['coefficients']}"
                    for i in range(len(shell['coefficients'])-1):
                        shell['coefficients'][i][-1]=f"0.0"
                shell['coefficients']=[ [ f"{round(float(y),6):.7f}" for y in x ] for x in shell['coefficients'] ]
            new_shells.append(shell)

        new_bas['elements'][atty]['electron_shells']=new_shells
    return new_bas

def make_psi4_block(basis_dict:dict, basis_name='custom') -> str:
    new_bs_str_base=get_str(basis_dict,'psi4')
    ind=4*' '
    basis_blocs=[]
    new_bs_str=[f"basis {basis_name.replace('-','_')} {{",
                ind+f"# Basis set obtained by augmenting"]+[
        ind+f"assign {basis_name}_custom", ""
    ]
    #for atty,b in parts.items():
    new_bs_str+=[ ind+f"[{basis_name}_custom]" ]
    new_bs_str+=[ ind+f"{x}" for x in new_bs_str_base.splitlines() if not x.strip()=='']
    new_bs_str+=[f"}}"]
    return new_bs_str

def produce_augmentation(basis_name, elements:List[str]):
    bs_dict=dict()
    for el in elements:
        if not requires_augmentation(basis_name):
            bs_dict.update(get_basis(basis_name, elements=elements))
        else:
            # First character indicates augmentation level (has to be [d,t,q]-aug-) as checked by requires_augmentation
            level = basis_name[0]
            aug_level=augmentation_map[level]
            base = basis_name.replace(f'{level}-aug-','aug-')
            # For H the aug-cc-pc basis is not available (since there is no core) need to fix that here
            if el.upper()=='H':
                if re.match(r'.*cc-pcv',basis_name.lower()):
                    base=base.lower().replace('cc-pcv','cc-pv')
            
            old_bs_dict = get_basis(base, elements=elements)
            bs_dict.update(geom_aug(old_bs_dict, aug_level, use_copy=True))

    # Proper formatting (and removal of redundant contractions)
    new_bs_dict=modify_basis(bs_dict)

    return new_bs_dict
#     def split_basis(basis_str):
#         lines=basis_str.splitlines()
#         spliter=r'[\*]{3,}'
#         import re
#         indices=[i for i,x in enumerate(lines) if re.match(spliter,x)]
#         basis_parts={}
# 
#         start=1
#         for i,idx in enumerate(indices):
#             if i>=start and (i-start)%1==0:
#                 header=lines[indices[i-1]+1]
#                 try:
#                     atty, z = header.split()
#                     int(z)
#                 except:
#                     print(f"Error parsing header line: {header}")
#                     continue
#                 block=lines[indices[i-1]+2:idx]
# 
#                 basis_parts[atty] = block
#         return basis_parts