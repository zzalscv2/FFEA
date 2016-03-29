import sys, os
import FFEA_surface

if len(sys.argv) != 2:
	sys.exit("Usage: python surf_load.py [INPUT FFEA surface (.surf)]")

fname = sys.argv[1]

surf = FFEA_surface.FFEA_surface(fname)
surf.print_details()
