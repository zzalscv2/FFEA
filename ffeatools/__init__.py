#from .FFEA_initialise import *
import enthusiasm as _enthusiasm
import os as _os
import sys as _sys

#add the modules folder (temporary hack until I can refactor)

_sys.path.append(_os.path.join(_os.path.dirname(_enthusiasm.__file__), 'modules'))
    
#from .FFEA_initialise import convert_pdb_to_emdb_map.py # need to make python wrapper first
#from .FFEA_initialise imp
    
from . import modules