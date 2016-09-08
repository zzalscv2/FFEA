#from .FFEA_initialise import *
import enthusiasm as _enthusiasm
import os as _os
import sys as _sys

#Many of these scripts expect that the modules folder is on the system or python path.
_sys.path.append(_os.path.join(_os.path.dirname(_enthusiasm.__file__), 'modules'))
    
#from . import modules
#from . import FFEA_initialise
#from . import FFEA_analysis

from .modules import __init__
from .FFEA_initialise import __init__
from .FFEA_analysis import __init__
