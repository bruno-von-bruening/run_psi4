
#!/usr/bin/env python
import re

def make_mp2_options():
    """ DF-MP2 settings as of Psi4"""
    return dict(
        mp2_type='\'df\'',
        freeze_core='True',
    )
def make_cc_options():
    """ DF-MP2 settings as of Psi4"""
    return dict(
        cc_type='\'df\'',
        freeze_core='True',
    )
def make_hf_options(open_spin:bool=False):
    """ DF-MP2 settings as of Psi4"""
    # https://psicode.org/psi4manual/1.9.x/scf.html
    return dict(
        reference='\'rhf\'' if not open_spin else '\'uhf\'',
        scf_type='\'direct\'',
    )
def make_dft_options(open_spin:bool=False):
    """ DF-MP2 settings as of Psi4"""
    return dict(
        reference='\'rks\'' if not open_spin else '\'uks\''
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

    cc  = True if match_method('cc',string)  else False
    mp2 = True if match_method('mp2',string)  else False
    hf = True if match_method('hf',string)  else False
    return {
        'mp2': mp2,
        'cc':  cc,
        'hf':  hf,
    }



if __name__=='__main__':
    raise Exception(f"Just a module at the moment")