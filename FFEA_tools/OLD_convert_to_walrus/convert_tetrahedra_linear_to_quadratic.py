import sys

def find_edge(n1, n2, edge_list):
	for edge in edge_list:
		if edge[0] == n1 and edge[1] == n2:
			return edge[2]

if len(sys.argv) != 3:
        sys.exit("Usage: python convert_tetrahedra_linear_to_quadratic [NODE FILE] [TOPOLOGY FILE]")

print "Running: convert_tetrahedra_linear_to_quadratic.py"

nodefile = open(sys.argv[1], "r")
topfile = open(sys.argv[2], "r")

# check that this is indeed a walrus node file
check = nodefile.readline().rstrip()
print "file type = " + check
if check != 'walrus node file':
	sys.exit("Error: walrus node file is expected for argument 1")

nodefile.readline() # skip "num_nodes" line
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

topfile.readline() # skip "num_elements" line
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

# create the edge list
edge_list = []
for el in elements:
	# sort the indices in ascending order
	e = el.split()
	e.sort()

	# get the node indices
	n1 = e[0]
	n2 = e[1]
	n3 = e[2]
	n4 = e[3]

	# there are 6 edges in a linear tetrahedron
	edge_list.append((n1, n2, -1))
	edge_list.append((n1, n3, -1))
	edge_list.append((n1, n4, -1))
	edge_list.append((n2, n3, -1))
	edge_list.append((n2, n4, -1))
	edge_list.append((n3, n4, -1))


# remove duplicate edges from list
edge_list = [[edge[0], edge[1], edge[2]] for edge in list(set(edge_list))]

for edge in edge_list:
	# set the index of the new midpoint node in the edge list
	edge[2] = str(len(nodes))

	# calc and add new midpoint node to the node list for the current edge
	node_a = nodes[int(edge[0])].split()
	node_b = nodes[int(edge[1])].split()
	new_node_x = (float(node_a[0]) + float(node_b[0])) * .5
	new_node_y = (float(node_a[1]) + float(node_b[1])) * .5
	new_node_z = (float(node_a[2]) + float(node_b[2])) * .5
	nodes.append(str(new_node_x) + " " + str(new_node_y) + " " + str(new_node_z) + "\n")

# construct new element topology
new_elements = []
for el in elements:
	# sort the indices in ascending order
	e = el.split()
	e.sort()

	# get the node indices
	n1 = e[0]
	n2 = e[1]
	n3 = e[2]
	n4 = e[3]

	new_el = n1 + " " + n2 + " " + n3 + " " + n4
	new_el += " " + find_edge(n1, n2, edge_list)
	new_el += " " + find_edge(n1, n3, edge_list)
	new_el += " " + find_edge(n1, n4, edge_list)
	new_el += " " + find_edge(n2, n3, edge_list)
	new_el += " " + find_edge(n2, n4, edge_list)
	new_el += " " + find_edge(n3, n4, edge_list)

	new_elements.append(new_el)


# write out new node list
nodefile = open(sys.argv[1], "w")
nodefile.write('walrus node file\n')
nodefile.write('num_nodes ' + str(len(nodes)) + "\n")
nodefile.write('num_surface_nodes ?\n')
nodefile.write('num_interior_nodes ?\n')
nodefile.write('surface nodes:\n')
nodefile.write('interior nodes:\n')
for node in nodes:
	nodefile.write(node)
nodefile.close()

# write out new topology
topfile = open(sys.argv[2], "w")
topfile.write('walrus topology file\n')
topfile.write('num_elements ' + str(len(elements)) + "\n")
topfile.write('num_surface_elements ?\n')
topfile.write('num_interior_elements ?\n')
topfile.write('surface elements:\n')
topfile.write('interior elements:\n')

for el in new_elements:
	topfile.write(el + "\n")
topfile.close()

print "Done. --> convert_tetrahedra_linear_to_quadratic.py"
