
#!/usr/bin/env python
import re

from enum import Enum
class method_types(str,Enum):
    MP = 'mp'
    CC = 'cc'
    HF = 'hf'
    DFT = 'dft'
    CI = 'ci'

def decorate_scf_type(scf_type:str):
    return "'{}'".format(scf_type.strip('\'\"'))


def post_hf_options(tag:method_types, scf_type:str):
    # https://psicode.org/psi4manual/master/autodoc_glossary_options_c.html#term-CC_TYPE-GLOBALS
    # https://psicode.org/psi4manual/1.4.0/autodir_options_c/globals__mp_type
    """ For MP and CC methods """
    if scf_type in ['cd','df']:
        return {
                f"{tag}_type": decorate_scf_type(scf_type),
                f"scf_type": decorate_scf_type(scf_type),
                }
    elif scf_type in ['conv']:
        return {
                f"{tag}_type": decorate_scf_type(scf_type),
                f"scf_type": decorate_scf_type('pk'),
        }
    else: raise ValueError(f"scf_type '{scf_type}' not recognized for {tag}!")

def make_mp_options(order:int, scf_type:str):
    """ DF-MP settings as of Psi4"""
    return post_hf_options(f"mp{order}", scf_type=scf_type)

def make_mp2_options(scf_type:str='default'):
    return make_mp_options(2,scf_type=scf_type)
def make_mp3_options(scf_type:str='default'):
    return make_mp_options(3,scf_type=scf_type)
def make_mp4_options(scf_type:str='default'):
    return make_mp_options(4,scf_type=scf_type)
def make_mp_options_gen(method:str):
    """ DF-MP settings as of Psi4"""
    if method.lower()=='mp2':
        func=make_mp2_options
    elif method.lower()=='mp3':
        func=make_mp3_options
    elif method.lower()=='mp4':
        func=make_mp4_options
    else: raise Exception(f"Unkown method {method}")
    return func

def make_cc_options(scf_type:str):
    """ DF-MP2 settings as of Psi4"""
    return post_hf_options("cc", scf_type=scf_type)
def make_ci_options(scf_type:str):
    """ CI settings as of Psi4"""
    assert scf_type in ['conv'], f"scf_type {scf_type} not recognized for CI!"
    return dict(
        ci_type=(decorate_scf_type(scf_type), "There should be a df option but it raises and error when I try to set"),
    )
def make_hf_options(scf_type:str): # https://psicode.org/psi4manual/1.4.0/autodir_options_c/globals__scf_type
    """ DF-MP2 settings as of Psi4"""
    # https://psicode.org/psi4manual/1.9.x/scf.html
    assert scf_type in ['direct','df','pk','mem_df','disk_df','out_of_core','cd','gtfock'], f"scf_type {scf_type} not recognized!"
    return dict(
        scf_type=decorate_scf_type(scf_type)
    )
def make_dft_options(scf_type:str):
    """ DF-MP2 settings as of Psi4"""
    return make_hf_options(scf_type)

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
