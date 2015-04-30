# Reads through the surface file to discover which nodes lie on the surface. Then calculates a
# mapping for the node indices such that the surface nodes all lie at the start of the node list.
# pyrar@leeds.ac.uk

import sys

def add_node(surface_node_list, a):
	for n in surface_node_list:
		if n == a:
			return
	surface_node_list.append(a)
	return


if len(sys.argv) != 4:
	sys.exit("Usage: python put_surface_nodes_first.py [NODE FILE] [TOPOLOGY FILE] [SURFACE FILE]")

print "Running: put_surface_nodes_first.py"

nodefile = open(sys.argv[1], "r")
topfile = open(sys.argv[2], "r")
surffile = open(sys.argv[3], "r")

# check that this is indeed a walrus node file
check = nodefile.readline().rstrip()
print "file type = " + check
if check != 'walrus node file':
	sys.exit("Error: walrus node file is expected for argument 1")

check = nodefile.readline().rstrip().split()
num_nodes = int(check[1])

nodefile.readline() # skip "num_surface_nodes" line
nodefile.readline() # skip "num_interior_nodes" line
nodefile.readline() # skip "surface nodes:" line

check = nodefile.readline().rstrip()
if check != 'interior nodes:':
	sys.exit("Error: Line 5 of node file should read 'interior nodes:'. Node file is either badly constructed, or has already been run through later scripts...")

# read in all the nodes
nodes = nodefile.readlines()

# check that this is indeed a walrus topology file
check = topfile.readline().rstrip()
print "file type = " + check
if check != 'walrus topology file':
	sys.exit("Error: walrus topology file is expected for argument 2")

check = topfile.readline().rstrip().split()
num_elements = int(check[1])

topfile.readline() # skip "num_surface_elements" line
topfile.readline() # skip "num_interior_elements" line
topfile.readline() # skip "surface elements:" line

check = topfile.readline().rstrip()
if check != 'interior elements:':
	sys.exit("Error: Line 5 of topology file should read 'interior elements:'. Topology file is either badly constructed, or has already been run through later scripts...")

# read in all the elements
elements = topfile.readlines()

nodefile.close()
topfile.close()


# check that this is indeed a walrus surface file
check = surffile.readline().rstrip()
print "file type = " + check
if check != 'walrus surface file':
	sys.exit("Error: walrus surface file is expected for argument 2")

check = surffile.readline().rstrip().split()
num_faces = int(check[1])

check = surffile.readline().rstrip()
if check != 'faces:':
	sys.exit("Error: Line 3 of surface file should read 'faces:'.")

# read in all the faces to a big list
faces = surffile.readlines()

surffile.close()

print "num nodes in node file = " + str(num_nodes)
print "num elements in topology file = " + str(num_elements)
print "num faces in surface file = " + str(num_faces)


# create an empty list to which the surface nodes can be added
surface_node_list = []

print "Creating list of surface nodes..."
for i in range(num_faces):
	a = faces[i].split()
	a[1] = int(a[1])
	a[2] = int(a[2])
	a[3] = int(a[3])
	add_node(surface_node_list, a[1])
	add_node(surface_node_list, a[2])
	add_node(surface_node_list, a[3])

num_surface_nodes = len(surface_node_list)
print "done. Found " + str(num_surface_nodes) + " surface nodes."

# get the list of non-surface nodes
print "Creating list of non-surface nodes..."
non_surface_node_list = range(num_nodes)

for sn in surface_node_list:
	non_surface_node_list.remove(sn);

num_non_surface_nodes = len(non_surface_node_list)
print "done. " + str(num_non_surface_nodes) + " non-surface nodes."

# create a list to hold the node swap mapping. Start with identity map.
print "Creating node index mapping..."
map = range(num_nodes)
lolfrus = range(num_nodes)

# Now change this mapping so that it maps surface nodes to a block at the start of the node list:
for i in range(num_surface_nodes):
	map[surface_node_list[i]] = i

# And the non surface nodes to a block at the end:
for i in range(num_non_surface_nodes):
	map[non_surface_node_list[i]] = i + num_surface_nodes

# Apply the mapping to the topology file
topfile = open(sys.argv[2], "w")
topfile.write('walrus topology file\n')
topfile.write('num_elements ' + str(num_elements) + '\n')
topfile.write('num_surface_elements ?\n')
topfile.write('num_interior_elements ?\n')
topfile.write('surface elements:\n')
topfile.write('interior elements:\n')
for el in elements:
	a = el.split()
	for j in range(10):
		a[j] = int(a[j])
		topfile.write(str(map[a[j]]) + " ");
	topfile.write("\n");
topfile.close()

# Apply the mapping to the surface file
surffile = open(sys.argv[3], "w")
surffile.write('walrus surface file\n')
surffile.write('num_surface_faces ' + str(num_faces) + '\n')
surffile.write('faces:\n')
for f in faces:
	a = f.split()
	a[1] = int(a[1])
	a[2] = int(a[2])
	a[3] = int(a[3])
	surffile.write(a[0] + " " + str(map[a[1]]) + " " + str(map[a[2]]) + " " + str(map[a[3]]) + "\n");
surffile.close()

# Apply the mapping to the node list (note that the sort step destroys the map structure so it can't be used after this)
paired = zip(map, nodes)
paired.sort()

nodefile = open(sys.argv[1], "w")
nodefile.write('walrus node file\n')
nodefile.write('num_nodes ' + str(num_nodes) + "\n");
nodefile.write('num_surface_nodes ' + str(num_surface_nodes) + "\n");
nodefile.write('num_interior_nodes ' + str(num_nodes - num_surface_nodes) + "\n");
nodefile.write('surface nodes:\n');
count = 0
for n in paired:
	nodefile.write(n[1])
	count += 1
	if count == num_surface_nodes:
		nodefile.write('interior nodes:\n');

nodefile.close()

print "Done. --> put_surface_nodes_first.py"
