import sys, os
import FFEA_node

if (len(sys.argv) != 6):
	sys.exit("Usage: python FFEA_node_rotate.py [INPUT .node fname] [OUTPUT .node fname] [x] [y] [z]")
	
# Get args
infname = sys.argv[1]
outfname = sys.argv[2]
R = [float(i) for i in sys.argv[3:]]
node = FFEA_node.FFEA_node(infname)
node.rotate(R)
node.write_to_file(outfname)