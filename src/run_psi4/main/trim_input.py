from . import *

def psi4_basis_set_escapes(basis_set):
    """ Psi4 escapes special charactes in basis set definitions:
    See: https://psicode.org/psi4manual/1.9.x/basissets.html"""
    escapes={
        r'\*':'s',
        r'\+':'p',
        r'[()\,]':'_',
    }
    to_escape=[]
    for escapable, replacement in escapes.items():
        the_search=re.findall(escapable, basis_set)
        to_escape+=[ (k, replacement) for k in the_search]
    
    if len(to_escape)>0:
        my_info(f"Found characters that psi4 will need to escape in basis set {basis_set}: {[x[0] for x in to_escape]}")
        for tag, rep in to_escape:
            basis_set=basis_set.replace(tag,rep)

    return basis_set


def check(input,object_type, default={}, ):
    """ Care about defaults """
    if isinstance(input, list): 
        assert all([isinstance(x,str) for x in func])

    elif isinstance(input, str):
        return input
    
    elif isinstance(input, type(None)):
        if not do_test:
            input=default['default_test']
            print(f"No {object_type} provided ressort to default: {input}")
        else:
            input=default['default']
    else:
        
        raise Exception()
    return input

def check_func_known(func):
    func_def={
        'default':['pbe0'],
        'default_test':['PBE0','PBE0-GRAC','WB97X-iptuning','MP2', 'CCSD']
    }
    return check(func,'functional', func_def)

def check_basis_known(basis_set):
    basis_def={
        'default':['sto-3g'],
        'default_test':['sto-3g']
    }
    return check(basis_set,'basis_set', basis_def)