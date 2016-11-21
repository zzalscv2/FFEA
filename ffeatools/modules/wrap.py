# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 09:39:41 2016

@author: py12rw
"""

import subprocess as _subprocess

def wrap_process(name, argv):
    args = [name]
    for arg in argv:
        args.append(arg)
    print("Trying args: "+str(args))
    _subprocess.call(args)
    return
    
def sanitize_tuple(args):
    sanitized = ()
    for arg in args:
        sanitized = sanitized + (str(arg),)
    return sanitized