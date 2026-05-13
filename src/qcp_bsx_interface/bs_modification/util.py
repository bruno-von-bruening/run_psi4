from qcp_bsx_interface.util.import_helper import *


def format_basis(basis_dict):
    import copy
    new_bas=copy.deepcopy(basis_dict)
    elements=new_bas['elements']
    psi4_formatted=lambda x: f"{round(float(x),5):.7f}"

    for atty,el_data in elements.items():
        vals=el_data['electron_shells']
        new_shells=[]

        for i,shell in enumerate(vals):

            if 'exponents' in shell:
                shell['exponents']=[ psi4_formatted(x) for x in shell['exponents'] ]

            if 'coefficients' in shell:
                ## This inirLLY DECOUPLE the last entry from the contracted basis set (not generally valid (e.g. core polarized)) 
                ## if len(shell['coefficients'])>1:
                ##     #assert len(shell['coefficients'])==2, f"Unexpected number of contractions in shell (expected 2 one for flexible diffuse, one for core): {len(shell['coefficients'])}: {shell['coefficients']}"
                ##     for i in range(len(shell['coefficients'])-1):
                ##         print(shell['coefficients'][i][-1])
                ##         shell['coefficients'][i][-1]=f"0.0"
                
                shell['coefficients']=[ [ psi4_formatted(y) for y in x ] for x in shell['coefficients'] ]
            new_shells.append(shell)

        new_bas['elements'][atty]['electron_shells']=new_shells
    return new_bas
formatted_basis_dict=Annotated[ Dict[str, Any], BeforeValidator(format_basis)] # for better readability of type hints

HIGHER_AUG_REGEX=r'[d,t,q]-aug'
AUGMENTATION_MAP = {
    'd': 1,
    't': 2,
    'q': 3,
}
KNOWN_PARTICLES=dict(
    conditional_core_polarization = [
        r'\-wc\drfnh',r'\-c\drfnh',
        r'\-wc\dr', r'\-c\dr',
    ], # integer tells which row # THE MINUS is essential (otherwise latter may match first!)
    augmentation=[HIGHER_AUG_REGEX,'aug'],
    dunning_stem=[r'cc-pv(wc|c)?[d,t,q,5,6]z'],
    ri=[r'\-ri']
)

@val_call
def basis_recognized(basis_name:str) -> str:
    # Check if the basis name matches the pattern for augmentation (should consist of known particles separated by dashes)
    remainder=basis_name.lower()
    for regex in KNOWN_PARTICLES.values():
        if not isinstance(regex, list):
            regex = [regex]

        for r in regex:
            if bool(re.search(r, remainder)):
                remainder = re.sub(r, '', remainder)
    # remainder should now only consist of dashes
    if not all(c == '-' for c in remainder):
        raise ValueError(f"Basis name {basis_name} has unkown particles: {remainder}")
    else:
        return basis_name
recognized_basis=Annotated[str, BeforeValidator(basis_recognized)]