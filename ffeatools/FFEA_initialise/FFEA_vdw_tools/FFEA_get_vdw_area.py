import sys
import FFEA_surface, FFEA_node
import FFEA_vdw

import argparse as _argparse

parser = _argparse.ArgumentParser(description="Get VDW Area")
parser.add_argument("vdwfname", action="store", help="Input .vdw file")
parser.add_argument("surffname", action="store", help="Input .surf file")
parser.add_argument("nodefname", action="store", help="Input .node file")

def FFEA_get_vdw_area(vdwfname, surffname, nodefname, print_only=False):
    # Open files
    vdw = FFEA_vdw.FFEA_vdw(vdwfname)
    surf = FFEA_surface.FFEA_surface(surffname)
    node = FFEA_node.FFEA_node(nodefname)
    
    areas = vdw.calc_active_areas(surf, node)
    VDW_index = []
    for i in range(-1,6,1):
        if print_only:
            print "VdW index " + str(i) + ": Area = " + str(areas[i + 1]) 
        else:
            VDW_index.append([i, areas[i+1]])
    if print_only:
        return
    return VDW_index
            

if sys.stdin.isatty():
    args = parser.parse_args()
    # Get args and build objects
    FFEA_get_vdw_area(args.vdwfname, args.surffname, args.nodefname, print_only=True)