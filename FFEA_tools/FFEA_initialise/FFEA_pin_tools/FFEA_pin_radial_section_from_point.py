import sys, os
import FFEA_node, FFEA_pin
import numpy as np
if len(sys.argv) != 8 and len(sys.argv) != 11:
	sys.exit("Usage: python FFEA_pin_radially.py [INPUT .node fname] [OUTPUT .pin fname] [lower radius (angstroms)] [upper radius (angstroms )] [Reference Point (x,y,z)] [Object centroid (x y z) {Optional}]")

# Get args
node_fname = sys.argv[1]
pin_fname = sys.argv[2]
low_radius = float(sys.argv[3])
up_radius = float(sys.argv[4])

# Get nodes and create pinned node object
nodes = FFEA_node.FFEA_node(node_fname)
pinned = FFEA_pin.FFEA_pin(pin_fname)

# If pinned node file exists, ask if we want to overwrite
if(os.path.exists(pin_fname)):
	cont = raw_input("File '" + pin_fname + "' already exists. Shall we overwrite? (y/n):")
	if(cont.lower()) == "y":
		pinned.reset()
	else:
		sys.exit("Please provide a new pinned filename.")

# Ref point
ref_point = np.array([float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])])

# Set nodal centroid
if len(sys.argv) == 11:
	centroid = [float(sys.argv[8]), float(sys.argv[9]), float(sys.argv[10])]
	nodes.set_centroid(centroid)

# Pin some nodes!
for i in range(nodes.num_nodes):
	distance = np.linalg.norm(nodes.pos[i] - ref_point)
	if distance <= up_radius and distance >= low_radius:
		pinned.add_node(i)

# Output
pinned.write_to_file(pin_fname)


