import sys, os
import FFEA_node

if len(sys.argv) != 2:
	sys.exit("Usage: python node_load.py [INPUT FFEA node (.node)]")

fname = sys.argv[1]

node = FFEA_node.FFEA_node(fname)
node.print_details()
