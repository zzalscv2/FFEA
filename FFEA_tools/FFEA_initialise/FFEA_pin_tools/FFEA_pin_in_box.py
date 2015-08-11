import sys
import FFEA_node, FFEA_pin
from numpy import linalg

if len(sys.argv) != 9:
	sys.exit("Usage: python FFEA_pin_in_box.py [INPUT .node fname] [OUTPUT .pin fname] [x limits (max min)] [y limits (max min)] [z limits (max min)]")

# Get args
node_fname = sys.argv[1]
pin_fname = sys.argv[2]
limits = [float(sys.argv[i + 3]) for i in range(6)]

# Get nodes and create pinned node object
nodes = FFEA_node.FFEA_node(node_fname)
pinned = FFEA_pin.FFEA_pin(pin_fname)
pinned.reset()

# Set nodal centroid
centroid = [0.0, 0.0, 0.0]
nodes.set_centroid(centroid)

# Pin some nodes!
for i in range(nodes.num_nodes):
	if nodes.pos[i][0] >= limits[0] and nodes.pos[i][0] <= limits[1] and nodes.pos[i][1] >= limits[2] and nodes.pos[i][1] <= limits[3] and nodes.pos[i][2] >= limits[4] and nodes.pos[i][2] <= limits[5]:
		pinned.add_node(i)

# Output
pinned.write_to_file(pin_fname)
