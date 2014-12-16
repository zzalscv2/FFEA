import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python DAT_to_walrus_stokes_file.py [NODE FILE] [OUTPUT STOKES FILE]")

input_node = open(sys.argv[1], "r")
output_stokes = open(sys.argv[2], "w")

#Reading nodes stuff
node_lines = input_node.readlines()
num_nodes = int(node_lines[1].split()[1])

#Get largest dimension
x_min = float("inf")
y_min = float("inf")
z_min = float("inf")
x_max = -1 * float("inf")
y_max = -1 * float("inf")
z_max = -1 * float("inf")

for i in range(6,6 + num_nodes):
	line = node_lines[i].split()
	if float(line[0]) < x_min:
		x_min = float(line[0])
	elif float(line[0]) > x_max:
		x_max = float(line[0])

	if float(line[0]) < y_min:
		y_min = float(line[0])
	elif float(line[0]) > y_max:
		y_max = float(line[0])

	if float(line[0]) < z_min:
		z_min = float(line[0])
	elif float(line[0]) > z_max:
		z_max = float(line[0])

dim = [(x_max - x_min) / 2.0, (y_max - y_min) / 2.0, (z_max - z_min) / 2.0]
order = [0,0,0]
if dim[0] >= dim[1] and dim[0] >= dim[2]:
	order[0] = 0
	if dim[1] > dim[2]:
		order[1] = 1
		order[2] = 2	
	else:
		order[1] = 2
		order[2] = 1

elif dim[1] > dim[0] and dim[1] > dim[2]:
	order[0] = 1
	if dim[0] > dim[2]:
		order[1] = 0
		order[2] = 2	
	else:
		order[1] = 2
		order[2] = 0
else:
	order[0] = 2
	if dim[0] > dim[1]:
		order[1] = 0
		order[2] = 1	
	else:
		order[1] = 1
		order[2] = 0

#Treat as big sphere
output_stokes.write("ffea stokes radii file\n")
s = ''
for i in xrange(num_nodes):
	s += str(dim[order[0]]/num_nodes) + "\n"

output_stokes.write("num_nodes " + str(num_nodes) + "\n")
output_stokes.write(s)
output_stokes.close()
