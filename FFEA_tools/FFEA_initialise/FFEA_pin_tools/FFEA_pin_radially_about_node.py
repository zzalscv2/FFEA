import FFEA_pin, FFEA_node
import sys

if len(sys.argv) != 5:
	sys.exit("Usage: python FFEA_pin_radially.py [INPUT .node] [OUTPUT .pin] [Node to pin about] [Radius]")

# Get args
node = FFEA_node.FFEA_node(sys.argv[1])
pin = FFEA_pin.FFEA_pin(sys.argv[2])
node_index = int(sys.argv[3])
radius = float(sys.argv[4])
pin.pin_radially(node, node.pos[node_index], radius)
pin.write_to_file(sys.argv[2])
