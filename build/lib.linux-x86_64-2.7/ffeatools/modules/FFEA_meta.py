# -*- coding: utf-8 -*-
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