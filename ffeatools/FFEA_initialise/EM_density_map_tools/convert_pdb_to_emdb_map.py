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
Created on Tue Aug 16 09:29:39 2016

@author: py12rw
"""

import wrap
import __builtin__

def convert_pdb_to_emdb_map(*argv):
    """
    This function is a wrapper for the convert_pdb_to_emdb_map executable.
    It can be used to convert atomistic PDB files to electron density maps. 
    ---
    Use it like this:
    convert_pdb_to_emdb_map(/path/to/input /path/to/output X_resolution Y_resolution Z_resolution atomic_radius)
    Where the atomic radius is given in angstroms.
    """

    # Log meta info    
    
    try:
        convert_pdb_to_emdb_map = {}
        convert_pdb_to_emdb_map["Input .PDB File"] = argv[0]
        convert_pdb_to_emdb_map["Output .map filename"] = argv[1]
        convert_pdb_to_emdb_map["X resolution"] = argv[2]
        convert_pdb_to_emdb_map["Y resolution"] = argv[3]
        convert_pdb_to_emdb_map["Z resolution"] = argv[4]
        convert_pdb_to_emdb_map["Atomic Radius (angstroms)"] = argv[5]
        __builtin__.meta_info.log["Convert DPB to EMDB map"] = convert_pdb_to_emdb_map
    except:
        print("Warning! failed to write log.")
        pass

    
    argv = wrap.sanitize_tuple(argv)
    wrap.wrap_process("./convert_pdb_to_emdb_map", argv)
    return
