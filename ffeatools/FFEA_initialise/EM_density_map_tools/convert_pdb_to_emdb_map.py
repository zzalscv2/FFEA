# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 09:29:39 2016

@author: py12rw
"""

import wrap

def convert_pdb_to_emdb_map(*argv):
    """
    This function is a wrapper for the convert_pdb_to_emdb_map executable.
    It can be used to convert atomistic PDB files to electron density maps. 
    ---
    Use it like this:
    convert_pdb_to_emdb_map(/path/to/input /path/to/output X_resolution Y_resolution Z_resolution atomic_radius)
    Where the atomic radius is given in angstroms.
    """
    argv = wrap.santitize_tuple(argv)
    wrap.wrap_process("./convert_pdb_to_emdb_map", argv)
    return