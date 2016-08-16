# Reads through the node file to discover which elements lie on the surface. Then calculates a
# mapping for the element indices such that the surface elements all lie at the start of the element list.
# pyrar@leeds.ac.uk

import sys

if len(sys.argv) != 5:
	sys.exit("Usage: python put_surface_elements_first.py [NODE FILE] [TOPOLOGY FILE] [SURFACE FILE] [MATERIALS FILE]")

print "Running: put_surface_elements_first.py"

nodefile = open(sys.argv[1], "r")
topfile = open(sys.argv[2], "r")
surffile = open(sys.argv[3], "r")
matfile = open(sys.argv[4], "r")

# check that this is indeed a walrus node file
check = nodefile.readline().rstrip()
print check
print "file type = " + check
if check != 'walrus node file' and check != 'ffea node file':
	sys.exit("Error: walrus/ffea node file is expected for argument 1")

nodefile.readline() # skip "num_nodes" line

# get the number of surface nodes
check = nodefile.readline().rstrip().split()
num_surface_nodes = int(check[1])

nodefile.close()

# check that this is indeed a walrus topology file
check = topfile.readline().rstrip()
print "file type = " + check
if check != 'walrus topology file' and check != 'ffea topology file':
	sys.exit("Error: walrus/ffea topology file is expected for argument 2")

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

# check that this is indeed a materials file
check = matfile.readline().rstrip()
print "file type = " + check
if check != 'walrus material params file' and 'ffea material params file':
	sys.exit("Error: walrus/ffea material params file is expected for argument 4")

check = matfile.readline().rstrip().split()
num_mat_elem = int(check[1])

# read in all the material elements
mat_elem = matfile.readlines()
print "len(mat_elem)", len(mat_elem)
#print mat_elem
matfile.close()


# check that this is indeed a walrus surface file
check = surffile.readline().rstrip()
print "file type = " + check
if check != 'walrus surface file' and check != 'ffea surface file':
	sys.exit("Error: walrus/ffea surface file is expected for argument 2")

check = surffile.readline().rstrip().split()
num_faces = int(check[1])

check = surffile.readline().rstrip()
if check != 'faces:':
	sys.exit("Error: Line 3 of surface file should read 'faces:'.")

# read in all the faces to a big list
faces = surffile.readlines()

surffile.close()


# create an empty list to which the surface elements can be added
surface_element_list = []
interior_element_list = []

print "Creating list of surface elements..."
for i in range(num_elements):
	a = elements[i].split()

	# if any of the elements nodes have indices below 'num_surface_nodes' then
	# they are classed as a "surface element" (an element containing at least 1 surface node)
	true = 0
	for node_index in a:
		if int(node_index) < num_surface_nodes:
			true = 1
			break

	if true == 1:
		surface_element_list.append(i)
	else:
		interior_element_list.append(i)


num_surface_elements = len(surface_element_list)
num_interior_elements = len(interior_element_list)
print "done. Found " + str(num_surface_elements) + " surface elements and " + str(num_interior_elements) + " interior elements."

# Now create a mapping of the form map[new_index] = old_index
map = surface_element_list + interior_element_list

# output the newly reordered topology file
topfile = open(sys.argv[2], "w")
topfile.write('ffea topology file\n')
topfile.write('num_elements ' + str(num_elements) + '\n')
topfile.write('num_surface_elements ' + str(num_surface_elements) + '\n')
topfile.write('num_interior_elements ' + str(num_interior_elements) + '\n')
topfile.write('surface elements:\n')
for i in range(num_elements):
	if i == num_surface_elements:
		topfile.write('interior elements:\n')
	topfile.write(elements[map[i]])	
if num_surface_elements == num_elements:
	topfile.write('interior elements:\n')

topfile.close()

# output the newly reordered materials file
matfile = open(sys.argv[4], "w")
matfile.write('ffea material params file\n')
matfile.write('num_elements ' + str(num_elements) + '\n')
print "len map = ", len(map)
print "len mat_elem = ", len(mat_elem)
for i in range(num_elements):
	matfile.write(mat_elem[map[i]])	
matfile.close()



# invert the mapping
inv_map = [[] for i in range(num_elements)]
for i in range(num_elements):
	inv_map[map[i]] = i

# remap the parent element index of surface faces and output the new surface file
surffile = open(sys.argv[3], "w")
surffile.write('ffea surface file\n')
surffile.write('num_surface_faces ' + str(num_faces) + '\n')
surffile.write('faces:\n')

for f in faces:
	f = f.split()
	f[0] = str(inv_map[int(f[0])])
	f = ' '.join(f)
	surffile.write(f + '\n')

surffile.close()

print "Done. --> put_surface_elements_first.py"
