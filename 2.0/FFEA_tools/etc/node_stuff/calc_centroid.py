import os, sys
from Vectors import *

if len(sys.argv) != 2:
	sys.exit("Usage: python calc_centroid.py [INPUT NODE FILE]")

node_lines = open(sys.argv[1], "r")

node_lines.readline() # ffea node file
num_nodes = int(node_lines.readline().split()[1])
num_surface_nodes = int(node_lines.readline().split()[1])
num_interior_nodes = int(node_lines.readline().split()[1])
node_lines.readline() # surface nodes:

centroid = vector3(0.0, 0.0, 0.0)
for i in range(num_surface_nodes):
	sline = node_lines.readline().split()
	centroid.x += float(sline[0])
	centroid.y += float(sline[1])
	centroid.z += float(sline[2])

node_lines.readline() # interior nodes:
for i in range(num_interior_nodes):
	sline = node_lines.readline().split()
	centroid.x += float(sline[0])
	centroid.y += float(sline[1])
	centroid.z += float(sline[2])

centroid.scale(1.0/num_nodes)
print centroid.x, centroid.y, centroid.z
