import pytest
from run_psi4.types.method import process_method

@pytest.mark.parametrize(
    "arg1,  exception",
    [
        # success cases: exception is None
        ('m06-hf', None),
        ('lrc-wpbe0-grac', None),
        ('hf', None),
        ('ccsd', None),
        # failure cases: expected is None
        (None, Exception),
        ('pbe0-asdf', Exception),
    ],
)
def test_process_method(arg1, exception):
    if exception:
        with pytest.raises(exception):
            process_method(arg1)
    else:
        process_method(arg1)
