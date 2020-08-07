# -*- coding: utf-8 -*-
# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

"""
Created on Tue Aug 16 09:39:41 2016

@author: py12rw
"""

import subprocess as _subprocess

def wrap_process(name, argv, bash=False):
    args = [name]
    for arg in argv:
        args.append(arg)
    print("Trying args: "+str(args))
    if bash:
        _subprocess.call(args, shell=True, executable="/bin/bash")
    else:
        _subprocess.call(args)
    return
    
def sanitize_tuple(args):
    sanitized = ()
    for arg in args:
        sanitized = sanitized + (str(arg),)
    return sanitized