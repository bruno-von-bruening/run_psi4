import pytest
from qcp_bsx_interface.bs_modification.process_basis import elements_require_core_polarization

@pytest.mark.parametrize("basis_name, elements, expected", [
    # wc3: apply core pol for elements in row 3+ (Z > 10), cutoff_z=10
    # S: Z=16 > 10 -> 'wc', H: Z=1 <= 10 -> None
    ('cc-pvdz',     ['O', 'H'], {'O':'cc-pvdz', 'H':'cc-pvdz'}),
    ('cc-pvdz-wc3r', ['S', 'H'], {'S': 'cc-pwcvdz', 'H': 'cc-pvdz'}),
    ('d-aug-cc-pvdz-wc3r', ['S', 'H'], {'S': 'd-aug-cc-pwcvdz', 'H': 'd-aug-cc-pvdz'}),
    # O: Z=8 <= 10 -> None, H: Z=1 <= 10 -> None
    ('cc-pvdz-wc3r', ['O', 'H'], {'O': 'cc-pvdz', 'H': 'cc-pvdz'}),
    ('cc-pvdz-wc2r', ['O', 'H'], {'O': 'cc-pwcvdz', 'H': 'cc-pvdz'}), # wc2r: apply core pol for elements in row 2+ (Z > 2), cutoff_z=2
    ('d-aug-cc-pvdz-wc2rfnh', ['S', 'Cl'], {'S': 'd-aug-cc-pwcvdz', 'Cl': 'd-aug-cc-pwcvdz'}), 
    ('d-aug-cc-pvdz-wc3rfnh', ['S', 'Cl'], {'S': 'd-aug-cc-pwcvdz', 'Cl': 'd-aug-cc-pvdz'}), 
    ('d-aug-cc-pvdz-wc4rfnh', ['S', 'Cl'], {'S': 'd-aug-cc-pvdz', 'Cl': 'd-aug-cc-pvdz'}), 
    # no core pol pattern -> None
])
def test_elements_require_core_polarization(basis_name, elements, expected):
    result = elements_require_core_polarization(basis_name, elements)
    assert result == expected

def test_unrecognized_basis_raises():
    with pytest.raises(ValueError):
        elements_require_core_polarization('cc-pvdz-wc3', ['S', 'H'])
