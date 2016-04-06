import sys, os
import FFEA_binding_sites

if len(sys.argv) != 2:
	sys.exit("Usage: python bsites_load.py [INPUT FFEA binding sites (.bsites)]")

fname = sys.argv[1]

bsites = FFEA_binding_sites.FFEA_binding_sites(fname)
bsites.print_details()
