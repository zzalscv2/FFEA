import sys, os
import numpy as np
import FFEA_node

if len(sys.argv) != 6:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT NODE FILE] [OUTPUT NODE FNAME] [DIRECTION (x/y/z)] [W.R.T (x/y/z)] [ANGLE (Degrees)]")

# Get args
innode_fname = sys.argv[1]
outnode_fname = sys.argv[2]

if sys.argv[3].lower() == "x":
	shear_dir = 0
elif sys.argv[3].lower() == "y":
	shear_dir = 1
elif sys.argv[3].lower() == "z":
	shear_dir = 2

if sys.argv[4].lower() == "x":
	wrt_dir = 0
elif sys.argv[4].lower() == "y":
	wrt_dir = 1
elif sys.argv[4].lower() == "z":
	wrt_dir = 2

angle = np.radians(float(sys.argv[5]))

nodes = FFEA_node.FFEA_node(innode_fname)

for i in range(len(nodes.pos)):
	nodes.pos[i][shear_dir] += nodes.pos[i][wrt_dir] / np.tan(angle)

nodes.write_nodes_to_file(outnode_fname)



