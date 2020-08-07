#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 05:23:33 2018

@author: rob
"""
try:
    import wrap
    import FFEA_script
    import FFEA_rod
except ImportError:
    from ffeatools import wrap
    from ffeatools import FFEA_script
    from ffeatools import FFEA_rod

def main():

    try:
        #wrap.wrap_process("ffeadev", ["symmetry_test_stretch_only.ffea"])
        wrap.wrap_process("../../../../src/ffea", ["symmetry_test_bend_only.ffea"])
        wrap.wrap_process("../../../../src/ffea", ["symmetry_test_stretch_only.ffea"])
    except OSError:
        raise
    
    bend_script = FFEA_script.FFEA_script("symmetry_test_bend_only.ffea")
    bend_rod = bend_script.rod[0]
    bend_rod.set_avg_energies()
    bend_analysis = FFEA_rod.anal_rod(bend_rod)
    bend_test_result = bend_analysis.do_bend_symmetry_test()
    
    stretch_script = FFEA_script.FFEA_script("symmetry_test_stretch_only.ffea")
    stretch_rod = stretch_script.rod[0]
    stretch_rod.set_avg_energies()
    stretch_analysis = FFEA_rod.anal_rod(stretch_rod)
    stretch_test_result = stretch_analysis.do_stretch_symmetry_test()
    
    if bend_test_result == False or stretch_test_result == False:
        raise SystemExit, 1

    raise SystemExit, 0



if __name__ == "__main__":
    main()