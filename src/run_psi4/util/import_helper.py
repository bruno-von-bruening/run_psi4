import re, os, sys, copy, time, shutil, time 
import inspect
import run_psi4.config.auto_options as p4ut
import numpy as np

from qcp_global_utils.encoding_and_conversion.encoding import nuclear_charge_to_element_symbol, element_symbol_to_nuclear_charge
from qcp_global_utils.pydantic.pydantic import pdtc_file, pdtc_directory, myBaseModel



# Data validation
from typing import Annotated, List, Union, Tuple,Dict
from pydantic import (
    BeforeValidator, 
    validate_call, 
    BaseModel, 
    Field
)
val_call=validate_call(config=dict(arbitrary_types_allowed=True),validate_return=True) # arbitraty types allowed to for non standard types

# Input Output
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
info=logger.info
warn=logger.warning
