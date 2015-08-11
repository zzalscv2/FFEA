import sys, os
import FFEA_node, FFEA_pin
from numpy import linalg

if len(sys.argv) != 4 and len(sys.argv) != 7:
	sys.exit("Usage: python FFEA_pin_radially.py [INPUT .node fname] [OUTPUT .pin fname] [radius (angstroms)] [Object centroid (x y z) {Optional}]")

# Get args
node_fname = sys.argv[1]
pin_fname = sys.argv[2]
radius = float(sys.argv[3])

# Get nodes and create pinned node object
nodes = FFEA_node.FFEA_node(node_fname)
pinned = FFEA_pin.FFEA_pin(pin_fname)
pinned.reset()

# If pinned node file exists, ask if we want to overwrite
if(os.path.exists(pin_fname)):
	cont = raw_input("File '" + pin_fname + "' already exists. Shall we overwrite? (y/n):")
	if(cont.lower()) == "y":
		pinned.reset()
	else:
		sys.exit("Please provide a new pinned filename.")

# Set nodal centroid
if len(sys.argv) == 7:
	centroid = [float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])]
	nodes.set_centroid(centroid)

# Pin some nodes!
for i in range(nodes.num_nodes):
	distance = linalg.norm(nodes.pos[i])
	if distance <= radius:
		pinned.add_node(i)

# Output
pinned.write_to_file(pin_fname)


