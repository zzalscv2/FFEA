import sys, os
from Vectors import vector3

if len(sys.argv) != 6:
	sys.exit("Usage python " + sys.argv[0] + " [INPUT .node file] [OUTPUT .node file] [x shift] [y shift] [z shift]")

# Get args
infile = sys.argv[1]
outfile = sys.argv[2]
translation = vector3(float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))

fin = open(infile, "r")
fout = open(outfile, "w")

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
