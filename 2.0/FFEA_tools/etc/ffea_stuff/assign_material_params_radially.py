import os, sys, math

if len(sys.argv) != 6:
	sys.exit("Usage: python assign_material_params_radially.py [INPUT .NODE FILE] [INPUT .MAT FILE] [OUTPUT .MAT FILE] [LOWER RADIUS LIMIT (angstroms)] [UPPER RADIUS LIMIT (angstroms)]")

if sys.argv[2] == sys.argv[3]:
	sys.exit("Specify different material files for input and output. Overwrite after if you want to.")
inputnode = open(sys.argv[1], "r")
inputmat = open(sys.argv[2], "r")
outputmat = open(sys.argv[3], "w")
radius_min = float(sys.argv[4])
radius_max = float(sys.argv[5])

nodelines = inputnode.readlines()
matlines = inputmat.readlines()
num_nodes = int(nodelines[1].split()[1])
num_surface_nodes = int(nodelines[2].split()[1])
num_interior_nodes = int(nodelines[3].split()[1])
nodes_to_change = []

for i in range(num_surface_nodes):
	line = nodelines[i + 5].split()
	radius = math.sqrt(float(line[0]) * float(line[0]) + float(line[1]) * float(line[1]))
	if radius > radius_min and radius < radius_max:
		nodes_to_change.append(i)
	
for i in range(num_interior_nodes):
	line = nodelines[6 + num_surface_nodes + i].split()
	radius = math.sqrt(float(line[0]) * float(line[0]) + float(line[1]) * float(line[1]))
	if radius > radius_min and radius < radius_max:
		nodes_to_change.append(i + num_surface_nodes)

#New parameters

print "Please enter:-\n"
density = raw_input("New density:")
shear_visc = raw_input("Shear viscosity:")
bulk_visc = raw_input("Bulk viscosity:")
shear_mod = raw_input("Shear modulus:")
bulk_mod = raw_input("Bulk modulus:")
dielec = raw_input("Dielectric constant:")

newmatline = str(float(density)) + " " + str(float(shear_visc)) + " " + str(float(bulk_visc)) + " " + str(float(shear_mod)) + " " + str(float(bulk_mod)) + " " + str(float(dielec)) + "\n"
index = -3
for matline in matlines:
	change = 0
	index = index + 1
	for node in nodes_to_change:
		if int(node) == index:
			change = 1
			break	
	if change == 1:
		outputmat.write(newmatline)
	else:
		outputmat.write(matline)
