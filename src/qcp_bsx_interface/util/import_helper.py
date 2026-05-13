from pydantic import validate_call; val_call=validate_call(config=dict(arbitrary_types_allowed=True, validate_return=True))


from typing import List, Dict, Tuple, Union, Annotated, Any
from pydantic import BeforeValidator, ValidationError
import  basis_set_exchange as bsx
import re
from basis_set_exchange  import write_formatted_basis_str as get_str, get_basis
from basis_set_exchange.manip import geometric_augmentation as geom_aug
from basis_set_exchange  import write_formatted_basis_str as get_str, get_basis

from qcp_global_utils.encoding_and_conversion.encoding import element_symbol_to_nuclear_charge, nuclear_charge_to_element_symbol
import copy, re