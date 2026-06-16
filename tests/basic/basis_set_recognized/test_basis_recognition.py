import subprocess as sp
import pytest
from contextlib import contextmanager
import os

description= \
"""
I recognize modification and raise and error if basis set is not recognized. Test here that standard basis sets will be
recognized by my code!
"""


@contextmanager
def jump(dir_name):
    origin = os.getcwd()   # 1. snapshot
    assert not '/' in dir_name.rstrip('/') # I do not want paths!
    
    location=os.path.dirname(__file__)
    path=os.path.join(location, dir_name)

    if os.path.isdir(path):
        p=sp.run(f"rm -r {path}/*", shell=True)
    else:
        os.mkdir(path)

    try:
        os.chdir(path)     # 2. move
        yield              # 3. hand control to `with` block
    finally:
        os.chdir(origin)   # 4. always return home


@pytest.mark.parametrize("bs", [
    'd-aug-cc-pvdz', 'd-aug-cc-pvcdz','d-aug-cc-pvwcdz', 'd-aug-pcseg-2',
    'def2-TZVPPD','def2-QZVP','def2-universal-jkfit',# ,'sto-3g',
    pytest.param('aug-cc-pcvdz',  marks=pytest.mark.xfail()),
])
def test_basis_recognized(bs):
    the_dir=f"{bs}_test"
    with jump(the_dir):
        cmd=f"run_psi4_exe.py --job_type polarizability_finite_field --method hf --geom dummy_h2 --basis {bs} --exc"
        print(f"Will run: {cmd}")
        p=sp.Popen(cmd, shell=True)
        p.communicate()
