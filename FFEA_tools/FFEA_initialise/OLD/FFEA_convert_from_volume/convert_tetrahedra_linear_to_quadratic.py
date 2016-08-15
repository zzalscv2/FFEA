import sys

def find_edge(n1, n2, edge_list):
	for edge in edge_list:
		if edge[0] == n1 and edge[1] == n2:
			return edge[2]

if len(sys.argv) != 4:
        sys.exit("Usage: python convert_tetrahedra_linear_to_quadratic [NODE FILE] [TOPOLOGY FILE] [STOKES FILE]")

print "Running: convert_tetrahedra_linear_to_quadratic.py"

nodefile = open(sys.argv[1], "r")
topfile = open(sys.argv[2], "r")
stokesfile = open(sys.argv[3], "r")

# check that this is indeed a walrus node file
check = nodefile.readline().rstrip()
print "file type = " + check
if check != 'walrus node file' and check != 'ffea node file':
	sys.exit("Error: walrus/ffea node file is expected for argument 1")

nodefile.readline() # skip "num_nodes" line
nodefile.readline() # skip "num_surface_nodes" line
nodefile.readline() # skip "num_interior_nodes" line
nodefile.readline() # skip "surface nodes:" line

check = nodefile.readline().rstrip()
if check != 'interior nodes:':
	sys.exit("Error: Line 5 of node file should read 'interior nodes:'. Node file is either badly constructed, or has already been run through later scripts...")

# read in all the nodes
nodes = nodefile.readlines()
nodefile.close()



# check that this is indeed a walrus stokes radii file
check = stokesfile.readline().rstrip()
print "file type = " + check
if check != 'walrus stokes radii file' and check != 'ffea stokes radii file':
	sys.exit("Error: walrus/ffea stokes radii file is expected for argument 3")

stokesfile.readline() # skip "num_nodes" line

# read in all the radii
radii = stokesfile.readlines()
stokesfile.close()



# check that this is indeed a walrus topology file
check = topfile.readline().rstrip()
print "file type = " + check
if check != 'walrus topology file' and check != 'ffea topology file':
	sys.exit("Error: walrus/ffea topology file is expected for argument 2")

topfile.readline() # skip "num_elements" line
topfile.readline() # skip "num_surface_elements" line
topfile.readline() # skip "num_interior_elements" line
topfile.readline() # skip "surface elements:" line

check = topfile.readline().rstrip()
if check != 'interior elements:':
	sys.exit("Error: Line 5 of topology file should read 'interior elements:'. Topology file is either badly constructed, or has already been run through later scripts...")

# read in all the elements
elements = topfile.readlines()
topfile.close()

# create the edge list
print "Creating edge list..."
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
print "Done."

print "Removing duplicates..."
# remove duplicate edges from list
edge_list = [[edge[0], edge[1], edge[2]] for edge in list(set(edge_list))]
print "Done."

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
	radii.append("0.0\n")

# construct new element topology
print "Building new elements..."
edge_dict = dict(((edge[0],edge[1]), edge[2]) for edge in edge_list)
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

#	new_el = n1 + " " + n2 + " " + n3 + " " + n4
#	new_el += " " + find_edge(n1, n2, edge_list)
#	new_el += " " + find_edge(n1, n3, edge_list)
#	new_el += " " + find_edge(n1, n4, edge_list)
#	new_el += " " + find_edge(n2, n3, edge_list)
#	new_el += " " + find_edge(n2, n4, edge_list)
#	new_el += " " + find_edge(n3, n4, edge_list)

	new_el = n1 + " " + n2 + " " + n3 + " " + n4
	new_el += " " + edge_dict[(n1, n2)]
	new_el += " " + edge_dict[(n1, n3)]
	new_el += " " + edge_dict[(n1, n4)]
	new_el += " " + edge_dict[(n2, n3)]
	new_el += " " + edge_dict[(n2, n4)]
	new_el += " " + edge_dict[(n3, n4)]

	new_elements.append(new_el)


print "Done."

# write out new node list
nodefile = open(sys.argv[1], "w")
nodefile.write('ffea node file\n')
nodefile.write('num_nodes ' + str(len(nodes)) + "\n")
nodefile.write('num_surface_nodes ?\n')
nodefile.write('num_interior_nodes ?\n')
nodefile.write('surface nodes:\n')
nodefile.write('interior nodes:\n')
for node in nodes:
	nodefile.write(node)
nodefile.close()

# write out new stokes radii file
stokesfile = open(sys.argv[3], "w")
stokesfile.write('ffea stokes radii file\n')
stokesfile.write('num_nodes ' + str(len(nodes)) + "\n")
for radius in radii:
	stokesfile.write(radius)
stokesfile.close()

# write out new topology
topfile = open(sys.argv[2], "w")
topfile.write('ffea topology file\n')
topfile.write('num_elements ' + str(len(elements)) + "\n")
topfile.write('num_surface_elements ?\n')
topfile.write('num_interior_elements ?\n')
topfile.write('surface elements:\n')
topfile.write('interior elements:\n')

for el in new_elements:
	topfile.write(el + "\n")
topfile.close()

print "Done. --> convert_tetrahedra_linear_to_quadratic.py"
