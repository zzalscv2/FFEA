import os, sys, math
from Vectors import vector3

if len(sys.argv) != 7:
	sys.exit("Usage: python vdw_radially.py [INPUT .NODE FILE] [INPUT .SURF FILE] [INPUT .VDW FILE] [OUTPUT .VDW FILE] [CUTOFF RADIUS (angstroms)] [VDW TYPE (-1 -> 6)]")

inputnode = open(sys.argv[1], "r")
inputsurf = open(sys.argv[2], "r")
inputvdw = open(sys.argv[3], "r")
outputvdw = open(sys.argv[4], "w")
radius = float(sys.argv[5])
new_vdw_type = sys.argv[6]

inlines = inputnode.readlines()
num_nodes = int(inlines[1].split()[1])
num_surface_nodes = int(inlines[2].split()[1])
num_interior_nodes = int(inlines[3].split()[1])
nodes = []

for i in range(num_surface_nodes):
	line = inlines[i + 5].split()
	nodes.append(line)
	
for i in range(num_interior_nodes):
	line = inlines[6 + num_surface_nodes + i].split()
	nodes.append(i + num_surface_nodes)


inlines = inputsurf.readlines()
if inlines[0] != "walrus surface file\n" and inlines[0] != "ffea surface file\n":
	sys.exit("Not a valid ffea surface file\n")

num_faces = int(inlines[1].split()[1])
faces = []
for i in range(num_faces):
	line = inlines[i + 3].split()
	faces.append(line)

inlines = inputvdw.readlines()
if int(inlines[1].split()[1]) != num_faces:
	sys.exit("Num faces in surface file and vdw file not equal")

vdw_types = []
for i in range(num_faces):
	line = inlines[i + 3]
	vdw_types.append(line)

outputvdw.write(inlines[0])
outputvdw.write(inlines[1])
outputvdw.write(inlines[2])

for i in range(num_faces):
	# node magnitude
	a_node = vector3(0.0, 0.0, 0.0)
	num_nodes_in_range = 0
	for j in faces[i]:
		a_node.x = float(nodes[int(j)][0])
		a_node.y = float(nodes[int(j)][1])
		a_node.z = float(nodes[int(j)][2])

		if a_node.mag() < radius:
			num_nodes_in_range += 1
	
	if num_nodes_in_range > 2:
		outputvdw.write(new_vdw_type + "\n")
	else:
		outputvdw.write(vdw_types[i])
	
