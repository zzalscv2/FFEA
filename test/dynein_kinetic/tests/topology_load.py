import sys, os
import FFEA_topology

if len(sys.argv) != 2:
	sys.exit("Usage: python topology_load.py [INPUT FFEA topology (.top)]")

fname = sys.argv[1]

top = FFEA_topology.FFEA_topology(fname)
top.print_details()
