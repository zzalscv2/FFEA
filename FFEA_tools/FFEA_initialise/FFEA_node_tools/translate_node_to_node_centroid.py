import sys, os
from math import *
from Vectors import *

if len(sys.argv) != 4:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT nodes_to_translate fname] [INPUT control_nodes_fname] [OUTPUT nodes_fname]")

translate_nodes_fname = sys.argv[1]
control_nodes_fname = sys.argv[2]
output_nodes_fname = sys.argv[3]

# Get centroids
centroid = [vector3(0.0, 0.0, 0.0)] * 2
fin = [open(translate_nodes_fname, "r"), open(control_nodes_fname, "r")]

for i in range(len(fin)):

	# Get header info
	fin[i].readline()
	num_nodes = int(fin[i].readline().split()[1])
	num_surface_nodes = int(fin[i].readline().split()[1])	
	num_interior_nodes = int(fin[i].readline().split()[1])

	# Surface nodes
	fin[i].readline()
	for j in range(num_surface_nodes):
		sline = fin[i].readline().split()
		centroid[i] += vector3(float(sline[0]), float(sline[1]), float(sline[2]))

	# Interior nodes
	fin[i].readline()
	for j in range(num_interior_nodes):
		sline = fin[i].readline().split()
		centroid[i] += vector3(float(sline[0]), float(sline[1]), float(sline[2]))
	
	fin[i].close()
	centroid[i].scale(1.0 / num_nodes)

# Get translation matrix
translation = centroid[1] - centroid[0]

# Move trans node to control nodes
fin = open(translate_nodes_fname, "r")
fout = open(output_nodes_fname, "w")

# Header stuff
fout.write(fin.readline())

# num_nodes
line = fin.readline()
num_nodes = int(line.split()[1])
fout.write(line)

# num_surface_nodes
line = fin.readline()
num_surface_nodes = int(line.split()[1])
fout.write(line)

# num_interior_nodes
line = fin.readline()
num_interior_nodes = int(line.split()[1])
fout.write(line)

# Write surface nodes
fout.write(fin.readline())
for i in range(num_surface_nodes):
	line = fin.readline()
	sline = line.split()
	node = vector3(float(sline[0]), float(sline[1]), float(sline[2]))
	node += translation
	fout.write("%6.3e %6.3e %6.3e\n" % (node.x, node.y, node.z))

# Write interior nodes
fout.write(fin.readline())
for i in range(num_interior_nodes):
	line = fin.readline()
	sline = line.split()
	node = vector3(float(sline[0]), float(sline[1]), float(sline[2]))
	node += translation
	fout.write("%6.3f %6.3f %6.3f\n" % (node.x, node.y, node.z))

fin.close()
fout.close()

