import os, sys, math

if len(sys.argv) != 4:
	sys.exit("Usage: python pin_radially.py [INPUT .NODE FILE] [OUTPUT .PIN FILE] [CUTOFF RADIUS (angstroms)]")

inputnode = open(sys.argv[1], "r")
outputpin = open(sys.argv[2], "w")
radius = float(sys.argv[3])

inlines = inputnode.readlines()
num_nodes = int(inlines[1].split()[1])
num_surface_nodes = int(inlines[2].split()[1])
num_interior_nodes = int(inlines[3].split()[1])
nodes_to_pin = []

for i in range(num_surface_nodes):
	line = inlines[i + 5].split()
	if math.sqrt(float(line[0]) * float(line[0]) + float(line[1]) * float(line[1])) < radius:
		nodes_to_pin.append(i)
	
for i in range(num_interior_nodes):
	line = inlines[6 + num_surface_nodes + i].split()
	if math.sqrt(float(line[0]) * float(line[0]) + float(line[1]) * float(line[1])) < radius:
		nodes_to_pin.append(i + num_surface_nodes)

outputpin.write("ffea pinned nodes file\n")
outputpin.write("num_pinned_nodes " + str(len(nodes_to_pin)) + "\n")
outputpin.write("pinned nodes:\n")
for node in nodes_to_pin:
	outputpin.write(str(node) + "\n")
