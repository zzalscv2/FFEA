import os, sys, math
import 
if len(sys.argv) != 9:
	sys.exit("Usage: python plot_stress.py [INPUT .NODE FILE] [INPUT STRESS .OUT FILE]")

inputnode = open(sys.argv[1], "r")
inputstress = open(sys.argv[2], "r")


