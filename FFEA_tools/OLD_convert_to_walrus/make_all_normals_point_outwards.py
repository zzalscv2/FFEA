# Calculates the normals for each face on the surface, and checks whether the normal points outwards or inwards
# (by checking against a point at the centre of the parent element). If it does not point outwards, then swap
# two indices of the face so that it does.
# pyrar@leeds.ac.uk

import sys

def in_or_out(x0, x1, x2, c):

	# calculate the normal of the face
	a = [x1[0] - x0[0], x1[1] - x0[1], x1[2] - x0[2]]
	b = [x2[0] - x0[0], x2[1] - x0[1], x2[2] - x0[2]]
	n = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

	# Calculate the projection of the normal with a vector from node 0 to the centroid of
	# the parent element. If it is negative then we know the normal is pointing outwards
	d = [c[0] - x0[0], c[1] - x0[1], c[2] - x0[2]]
	n_dot_d = n[0] * d[0] + n[1] * d[1] + n[2] * d[2]
	if n_dot_d < 0:
		return 1
	else:
		return 0


if len(sys.argv) != 4:
	sys.exit("Usage: python make_all_normals_point_outwards.py [NODE FILE] [TOPOLOGY FILE] [SURFACE FILE]")

print "Running: make_all_normals_point_outwards.py"

nodefile = open(sys.argv[1], "r")
topfile = open(sys.argv[2], "r")
surffile = open(sys.argv[3], "r")

# first check that this is indeed a walrus node file
check = nodefile.readline().rstrip()
print "file type = " + check
if check != 'walrus node file':
	sys.exit("Error: walrus node file is expected for argument 1")

# get total number of nodes
check = nodefile.readline().rstrip().split()
num_nodes = int(check[1])

# get number of surface nodes
check = nodefile.readline().rstrip().split()
num_surface_nodes = int(check[1])

# get number of interior nodes
check = nodefile.readline().rstrip().split()
num_interior_nodes = int(check[1])

# check for 'surface nodes:' line
check = nodefile.readline().rstrip()
if check != 'surface nodes:':
        sys.exit("Error: Line 5 of node file should read 'surface nodes:")

# put all the node positions into a big list
node = []
for i in range(num_surface_nodes):
	a = nodefile.readline().split()
	node.append([float(a[0]), float(a[1]), float(a[2])])

# check for 'interior nodes:' line
check = nodefile.readline().rstrip()
if check != 'interior nodes:':
        sys.exit("Error: expecting 'interior nodes:")

for i in range(num_interior_nodes):
	a = nodefile.readline().split()
	node.append([float(a[0]), float(a[1]), float(a[2])])


# check that we have a walrus topology file
check = topfile.readline().rstrip()
print "file type = " + check
if check != 'walrus topology file':
	sys.exit("Error: walrus topology file is expected for argument 2")

# get total number of elements
check = topfile.readline().rstrip().split()
num_elements = int(check[1])

# get number of surface elements
check = topfile.readline().rstrip().split()
num_surface_elements = int(check[1])

# get number of interior elements
check = topfile.readline().rstrip().split()
num_interior_elements = int(check[1])

# check for 'surface elements:' line
check = topfile.readline().rstrip()
if check != 'surface elements:':
        sys.exit("Error: Line 4 of topology file should read 'surface elements:")

# put all the element node indices into a big list
elem = []
for i in range(num_surface_elements):
	a = topfile.readline().split()
	elem.append([int(a[0]), int(a[1]), int(a[2]), int(a[3])])

# check for 'interior elements:' line
check = topfile.readline().rstrip()
if check != 'interior elements:':
        sys.exit("Error: expecting 'interior elements:")

for i in range(num_interior_elements):
	a = topfile.readline().split()
	elem.append([int(a[0]), int(a[1]), int(a[2]), int(a[3])])


# check that we have a walrus surface file
check = surffile.readline().rstrip()
print "file type = " + check
if check != 'walrus surface file':
	sys.exit("Error: walrus surface file is expected for argument 3")

# get total number of faces
check = surffile.readline().rstrip().split()
num_faces = int(check[1])

# check for 'faces:' line
check = surffile.readline().rstrip()
if check != 'faces:':
        sys.exit("Error: Line 3 of surface file should read 'faces:")

# put all faces in a big list
surface = []
for i in range(num_faces):
	a = surffile.readline().split()
	surface.append([int(a[0]), int(a[1]), int(a[2]), int(a[3])])

print "num nodes in node file = " + str(num_nodes)
print "num elements in topology file = " + str(num_elements)
print "num faces in surface file = " + str(num_faces)

nodefile.close()
topfile.close()
surffile.close()

print "Checking normals for all faces..."
outfile = open(sys.argv[3], "w")
outfile.write('walrus surface file\n')
outfile.write('num_surface_faces ' + str(num_faces) + '\n')
outfile.write('faces:\n')
for f in surface:
	el = elem[f[0]]
	n0 = node[el[0]]
	n1 = node[el[1]]
	n2 = node[el[2]]
	n3 = node[el[3]]
	c = [(n0[0] + n1[0] + n2[0] + n3[0])/4, (n0[1] + n1[1] + n2[1] + n3[1])/4, (n0[2] + n1[2] + n2[2] + n3[2])/4]
	if in_or_out(node[f[1]], node[f[2]], node[f[3]], c) == 1:
		outfile.write(str(f[0]) + " " + str(f[1]) + " " + str(f[2]) + " " + str(f[3]) + "\n")
	else:
		outfile.write(str(f[0]) + " " + str(f[1]) + " " + str(f[3]) + " " + str(f[2]) + "\n")

outfile.close()

print "Done. --> make_all_normals_point_outwards.py"
