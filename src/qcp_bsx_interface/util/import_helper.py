from pydantic import validate_call; val_call=validate_call(config=dict(arbitrary_types_allowed=True, validate_return=True))


from typing import List
import  basis_set_exchange as bsx
import re