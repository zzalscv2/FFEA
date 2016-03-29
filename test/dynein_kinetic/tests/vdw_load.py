import sys, os
import FFEA_vdw

if len(sys.argv) != 2:
	sys.exit("Usage: python vdw_load.py [INPUT FFEA vdw (.vdw)]")

fname = sys.argv[1]

vdw = FFEA_vdw.FFEA_vdw(fname)
vdw.print_details()
