import sys, os
import FFEA_pin

if len(sys.argv) != 2:
	sys.exit("Usage: python pin_load.py [INPUT FFEA pin (.pin)]")

fname = sys.argv[1]

pin = FFEA_pin.FFEA_pin(fname)
pin.print_details()
