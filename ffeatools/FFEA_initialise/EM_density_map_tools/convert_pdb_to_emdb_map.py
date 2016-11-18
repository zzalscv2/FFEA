# -*- coding: utf-8 -*-
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
