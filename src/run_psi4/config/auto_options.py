
#!/usr/bin/env python
import re

from enum import Enum
class method_types(str,Enum):
    MP = 'mp'
    CC = 'cc'
    HF = 'hf'
    DFT = 'dft'
    CI = 'ci'

def make_mp_options(use_df:bool=True):
    """ DF-MP settings as of Psi4"""
    return dict(
        mp_type="'df'" if use_df else "'conv'",
    )
def make_cc_options(use_df:bool=True):
    """ DF-MP2 settings as of Psi4"""
    return dict(
        cc_type="'df'" if use_df else "'conv'",
    )
def make_hf_options(use_df:bool=True):
    """ DF-MP2 settings as of Psi4"""
    # https://psicode.org/psi4manual/1.9.x/scf.html
    return dict(
        scf_type="'df'" if use_df else "'pk'",
    )
def make_dft_options(open_spin:bool=False, use_df:bool=True):
    """ DF-MP2 settings as of Psi4"""
    return dict(
        scf_type="'df'" if use_df else "'pk'",
    )
def make_ci_options():
    """ CI settings as of Psi4"""
    return dict(
        ci_type=("'conv'", "There should be a df option but it raises and error when I try to set"),
    )

def embed_value(x):
    """ In case its a string make it an explicit string (for printout) """
    if isinstance(x, str):
        return '\''+x+'\''
    else:
        return x
    
def get_method_flags(string):
    def match_method(tag, target):
        return bool(re.search(tag.upper(),target.upper()))

    matched=[]
    for method in method_types:
        if method == method_types.DFT: # DFT does not hat a leading tag
            continue
        elif match_method(method.value,string):
            matched+=[method.value]
    if len(matched)==0:
        method=method_types.DFT
    else:
        assert len(matched)<2, f"Multiple method tags matched in string '{string}': {matched}"
        method = matched[0]

    return method




if __name__=='__main__':
    raise Exception(f"Just a module at the moment")