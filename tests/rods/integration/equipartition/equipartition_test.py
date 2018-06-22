#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 05:23:33 2018

@author: rob
"""

# allow old-style pythonpath or new module imports
try:
    import wrap
    import FFEA_script
    import FFEA_rod
except ImportError:
    from ffeatools import wrap
    from ffeatools import FFEA_script
    from ffeatools import FFEA_rod
    
import numpy as np

def main():

    try:
        wrap.wrap_process("../../../../src/ffea", ["realistic.ffea"])
    except OSError:
        raise
    
    script = FFEA_script.FFEA_script("realistic.ffea")
    
    temp = 300 #kelvin
    analytical_kbT = temp*1.38064852 * 10**-23
    
    bendy_rod = script.rod[0]
    stretchy_rod = script.rod[1]
    twisty_rod = script.rod[2]
    
    #Bendy rod
    bendy_analysis = FFEA_rod.anal_rod(bendy_rod)
    bendy_analysis.get_equipartition()
    bendy_rod_avg_energy = np.average(bendy_analysis.bending_energy)
    print("Bendy energy = "+str(bendy_rod_avg_energy))
    print("Equipartition energy = "+str(0.5*analytical_kbT))
    bend_equal = FFEA_rod.rod_math.approximately_equal(0.5*analytical_kbT, bendy_rod_avg_energy, 0.08)
    print("Bendy status: "+str(bend_equal))
    print("---")
    
    #Stretchy rod
    stretchy_analysis = FFEA_rod.anal_rod(stretchy_rod)
    stretchy_analysis.get_average_quantities()
    stretchy_rod_avg_energy = stretchy_analysis.average_extension_sq*0.5*stretchy_analysis.get_constant_parameter(0)
    stretchy_time_avg_energy = np.average(stretchy_rod_avg_energy)
    print("Stretchy energy = "+str(stretchy_time_avg_energy))
    print("Equipartition energy = "+str(0.5*analytical_kbT))
    stretch_equal = FFEA_rod.rod_math.approximately_equal(0.5*analytical_kbT, stretchy_time_avg_energy, 0.08)
    print("Stretchy status: "+str(stretch_equal))
    print("---")
    
    #Twisty rod
    twisty_analysis = FFEA_rod.anal_rod(twisty_rod)
    twisty_analysis.get_equipartition()
    twisty_rod_avg_energy = np.average(twisty_analysis.twist_energy[:100])
    print("Tiwsty energy = "+str(twisty_rod_avg_energy))
    print("Equipartition energy = "+str(0.5*analytical_kbT))
    twist_equal = FFEA_rod.rod_math.approximately_equal(0.5*analytical_kbT, twisty_rod_avg_energy, 0.08)
    print("Twisty status: "+str(twist_equal))
    print("---")
    
    if twist_equal and bend_equal and stretch_equal:
        raise SystemExit, 0
    else:
        raise SystemExit, 1

if __name__ == "__main__":
    main()