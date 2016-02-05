import sys, os
import numpy as np
import FFEA_node, FFEA_kinetic_map

if len(sys.argv) != 4:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT .node fname] [OUTPUT .node fname] [INPUT ffea .map fname]\n")

# Get args
innode = sys.argv[1]
outnode = sys.argv[2]
inmap = sys.argv[3]

# Get nodes
input_nodes = FFEA_node.FFEA_node(innode)

# Get map
kinetic_map = FFEA_kinetic_map.FFEA_kinetic_map(inmap)

# Apply matrix!
output_nodes = kinetic_map.apply_sparse(input_nodes)

# Print to file

# Finding these values would take to long, so set all to surface as this is just a test script
num_nodes = len(output_nodes)
num_surface_nodes = len(output_nodes)
num_interior_nodes = 0

fout = open(outnode, "w")
fout.write("ffea node file\n")
fout.write("num_nodes %d\n" % (num_nodes))
fout.write("num_surface_nodes %d\n" % (num_surface_nodes))
fout.write("num_interior_nodes %d\n" % (num_interior_nodes))

# Surface nodes
fout.write("surface nodes:\n")
for i in range(num_surface_nodes):
	fout.write("%6.3f %6.3f %6.3f\n" % (output_nodes[i][0], output_nodes[i][1], output_nodes[i][2]))

# Interior nodes
fout.write("interior nodes:\n")
for i in range(num_surface_nodes, num_nodes, 1):
	fout.write("%6.3f %6.3f %6.3f\n" % (output_nodes[i][0], output_nodes[i][1], output_nodes[i][2]))

fout.close()
