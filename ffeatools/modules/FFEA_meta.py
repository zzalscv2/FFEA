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
Created on Fri Sep 23 10:03:43 2016

@author: py12rw
"""

import __builtin__

import json as _json
class FFEA_meta:
    def __init__(self):
        self.log = {}
        return
        
    def load_mega(self, fname):
        self.log = _json.loads(open(fname).read())

    def load_script(self, script):
        self.log["root_path"] = script.params.trajectory_out_fname.split("_trajectory.out")[0]

    def dump_log(self):
        try:
            out_filename = self.log["root_path"] + ".meta.json"
        except KeyError:
            out_filename = "FFEA_meta.json"
        with open(out_filename, 'w') as outfile:
            _json.dump(self.log, outfile)
            print("Dumped meta info to "+out_filename)
        import os
        
__builtin__.meta_info = FFEA_meta() # dirty hack for a good cause