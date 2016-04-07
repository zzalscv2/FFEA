import sys, os
import FFEA_stokes

if len(sys.argv) != 2:
	sys.exit("Usage: python stokes_load.py [INPUT FFEA stokes (.stokes)]")

fname = sys.argv[1]

stokes = FFEA_stokes.FFEA_stokes(fname)
stokes.print_details()
