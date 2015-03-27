import sys, os
from Vectors import vector3

if len(sys.argv) != 2:
	sys.exit("Usage python " + sys.argv[0] + " [INPUT .node file]")

# Get args
infile = sys.argv[1]
fin = open(infile, "r")

centroid = vector3(0.0, 0.0, 0.0)

# Header stuff
fin.readline()

# num_nodes
line = fin.readline()
num_nodes = int(line.split()[1])

# num_surface_nodes
line = fin.readline()
num_surface_nodes = int(line.split()[1])

# num_interior_nodes
line = fin.readline()
num_interior_nodes = int(line.split()[1])

# Get surface nodes
fin.readline()
for i in range(num_surface_nodes):
	line = fin.readline()
	sline = line.split()
	centroid += vector3(float(sline[0]), float(sline[1]), float(sline[2]))

# Get interior nodes
fin.readline()
for i in range(num_interior_nodes):
	line = fin.readline()
	sline = line.split()
	centroid += vector3(float(sline[0]), float(sline[1]), float(sline[2]))

fin.close()

centroid.scale(1.0 / num_nodes)
print "Centroid = (%6.3e, %6.3e, %6.3e)\n" % (centroid.x, centroid.y, centroid.z)
