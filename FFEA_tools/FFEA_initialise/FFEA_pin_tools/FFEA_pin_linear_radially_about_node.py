import FFEA_pin, FFEA_node, FFEA_topology
import sys

if len(sys.argv) != 6:
	sys.exit("Usage: python FFEA_pin_radially.py [INPUT .node] [INPUT .top] [OUTPUT .pin] [Node to pin about] [Radius]")

# Get args
node = FFEA_node.FFEA_node(sys.argv[1])
top = FFEA_topology.FFEA_topology(sys.argv[2])
pin = FFEA_pin.FFEA_pin(sys.argv[3])
node_index = int(sys.argv[4])
radius = float(sys.argv[5])
pin.pin_linear_radially(node, top, node.pos[node_index], radius)
pin.write_to_file(sys.argv[3])
