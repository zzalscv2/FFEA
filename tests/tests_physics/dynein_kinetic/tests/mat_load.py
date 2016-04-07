import sys, os
import FFEA_material

if len(sys.argv) != 2:
	sys.exit("Usage: python mat_load.py [INPUT FFEA material (.mat)]")

fname = sys.argv[1]

mat = FFEA_material.FFEA_material(fname)
mat.print_details()
