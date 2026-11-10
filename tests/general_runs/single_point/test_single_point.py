
from run_psi4.types.job import psi4_job
from run_psi4.main.main import psi4_run
import yaml, os, shutil
import pytest

try: # Only for python>=3.11 
    from contextlib import chdir
except Exception as ex:
    from contextlib import contextmanager
    @contextmanager
    def chdir(path):
        old_dir = os.getcwd()
        os.chdir(path)
        try:
            yield
        finally:
            os.chdir(old_dir)


THIS_DIR=os.path.dirname(__file__)
make_loc=lambda x:os.path.join(THIS_DIR,x)
@pytest.mark.parametrize(
    "conf_file",
    [
        make_loc('config_pbe0_df.yaml'),
        make_loc('config_MP2_df.yaml'),
        pytest.param(make_loc('config_MP2_conv.yaml'), marks=pytest.mark.xfail(reason="Conventional MP2 cannot evaluate density!")),
    ]
)
def test_run_single_point(conf_file):
    assert os.path.isfile(conf_file), f"Not a file: {conf_file}"

    new_dir=os.path.join( THIS_DIR, os.path.basename(conf_file).split('.')[0] )
    if not os.path.isdir(new_dir): os.mkdir(new_dir)

    conf_abs_path=os.path.realpath(conf_file)
    with chdir(new_dir):
        loc_conf_file=os.path.basename(conf_abs_path)
        shutil.copy(conf_abs_path,loc_conf_file)

        with open(loc_conf_file,'r') as rd:
            conf=yaml.safe_load(rd)
        conf.update(execute=True)

        job=psi4_job(**conf)

        psi4_run(job)
