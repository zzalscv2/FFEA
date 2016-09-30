# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 17:12:39 2016

@author: py12rw
"""

import sys

try:
    from ffeatools.modules import * # python package
except:
    try:
        import FFEA_trajectory
    except ImportError:
        print("Failure to import FFEA_trajectory")
        sys.exit(1) # failure to import

try:
    test_load_traj = FFEA_trajectory.FFEA_trajectory("unit_test_traj.ftj")
    sys.exit(0)
except IOError:
    print("Couldn't find trajectroy file that's supposed to be packed in with these tests. Probably a CMake issue!")
    sys.exit(1)
except Exception, e:
    print(e)
    sys.exit(1)