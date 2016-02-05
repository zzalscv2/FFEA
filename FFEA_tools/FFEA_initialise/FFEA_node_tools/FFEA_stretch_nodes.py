import sys, os
import FFEA_node
import numpy as np

if len(sys.argv) != 5:
	sys.exit("python FFEA_stretch_nodes.py [INPUT .node fname] [OUTPUT .node fname] [direction (x, y, z)] [Factor]")

# Get args
infname = sys.argv[1]
outfname = sys.argv[2]

if sys.argv[3].lower() == "x":
	direction = 0
elif sys.argv[3].lower() == "y":
	direction = 1
elif sys.argv[3].lower() == "z":
	direction = 2
else:
	sys.exit("Error. Direction not recognised.")

factor = float(sys.argv[4])

# Build nodes
node = FFEA_node.FFEA_node(infname)

# Get centroid
centroid = node.calc_centroid()

# All distances from centroid should be increased by a factor of 'factor' in the direction we require only
for p in node.pos:
	p[direction] = centroid[direction] + factor * (p[direction] - centroid[direction])

node.write_to_file(outfname)

