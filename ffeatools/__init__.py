#from .FFEA_initialise import *
from .modules import __init__
from .FFEA_initialise import __init__
from .FFEA_analysis import __init__

import enthusiasm as _enthusiasm
import os as _os
import sys as _sys
import signal as _signal
import __builtin__
import atexit as _atexit

#Many of these scripts expect that the modules folder is on the system or python path.
_sys.path.append(_os.path.join(_os.path.dirname(_enthusiasm.__file__), 'modules'))
    
#from . import modules
#from . import FFEA_initialise
#from . import FFEA_analysis

from .modules import FFEA_meta

def _save_meta(*argv):
    __builtin__.meta_info.dump_log()
    _os._exit(0)

for _sig in (_signal.SIGABRT, _signal.SIGINT, _signal.SIGTERM, _signal.SIGQUIT):
    _signal.signal(_sig, _save_meta)
    
_atexit.register(_save_meta)