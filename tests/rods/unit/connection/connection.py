#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 05:23:33 2018

@author: rob
"""

# allow old-style pythonpath or new module imports
import subprocess

def main():

    return_value = subprocess.call(["../../../../src/ffea", "connection.ffeatest"])
    raise SystemExit, return_value

if __name__ == "__main__":
    main()
