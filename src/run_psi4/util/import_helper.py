import re, os, sys, copy, time, shutil, time 
import inspect
import run_psi4.config.auto_options as p4ut
import numpy as np

from qcp_global_utils.encoding_and_conversion.encoding import nuclear_charge_to_element_symbol, element_symbol_to_nuclear_charge
from qcp_global_utils.pydantic.pydantic import pdtc_file, pdtc_directory, myBaseModel


from qcp_global_utils.parser.parser import controller_model

# Data validation
from typing import Annotated, List, Union, Tuple,Dict, Literal
from enum import Enum
from pydantic import (
    BeforeValidator, 
    validate_call, 
    BaseModel, 
    Field
)
val_call=validate_call(config=dict(arbitrary_types_allowed=True),validate_return=True) # arbitraty types allowed to for non standard types

# Input Output
import logging

# define a VERBOSE level between DEBUG(10) and INFO(20)
VERBOSE_LEVEL_NUM = 15
logging.addLevelName(VERBOSE_LEVEL_NUM, "VERBOSE")


def _verbose(self, message, *args, **kws):
    if self.isEnabledFor(VERBOSE_LEVEL_NUM):
        self._log(VERBOSE_LEVEL_NUM, message, args, **kws)


logging.Logger.verbose = _verbose


class ConditionalFormatter(logging.Formatter):
    """
    Use a compact format by default, and include file:line when level is VERBOSE (or DEBUG).
    """
    def __init__(self, fmt_normal, fmt_verbose, datefmt=None):
        super().__init__(fmt=fmt_normal, datefmt=datefmt)
        self.fmt_normal = fmt_normal
        self.fmt_verbose = fmt_verbose

    def format(self, record):
        # treat VERBOSE and DEBUG as verbose for formatting purposes
        if record.levelno == VERBOSE_LEVEL_NUM or record.levelno <= logging.DEBUG:
            self._style._fmt = self.fmt_verbose
        else:
            self._style._fmt = self.fmt_normal
        return super().format(record)


# create handler with conditional formatting
handler = logging.StreamHandler()
handler.setFormatter(
    ConditionalFormatter(
        fmt_normal="%(levelname)s: %(message)s",
        fmt_verbose="%(levelname)s: %(message)s (%(filename)s:%(lineno)d)",
    )
)

# configure module logger
logger = logging.getLogger(__name__)
logger.handlers = []  # avoid duplicate handlers if re-imported
logger.addHandler(handler)
logger.setLevel(logging.INFO)  # default; set to VERBOSE_LEVEL_NUM or logging.DEBUG for verbose output

info = logger.info
def warn(message, *args, **kwargs):
    logger.warning(message, *args, **kwargs)
    # print(f"WARNING: {message}")
