# Counts the number of occurences of triangular faces in the given topology file.
# Those faces which appear only once must lie on the surface. Surface faces are output
# to the specified surface file (preferred naming convention is "*.surf")
# pyrar@leeds.ac.uk

import sys

# First checks if the given face already exists. If so,
# it increments the counter for that face by 1. Otherwise,
# a new face is appended to the list with its counter set
# to 1.
def add_face(fl, n1, n2, n3, element):
	order = [n1, n2, n3]
	order.sort()
	n1 = order[0]
	n2 = order[1]
	n3 = order[2]
	for f in fl:
		if n1 == f[0] and n2 == f[1] and n3 == f[2]:
			f[3] += 1
			return
	fl.append([n1, n2, n3, 1, element])
	return


if len(sys.argv) != 3:
	sys.exit("Usage: python extract_surface_from_topology.py [INPUT TOPOLOGY FILE] [OUTPUT SURFACE FILE]")

print "Running: extract_surface_from_topology.py"

topfile = open(sys.argv[1], "r")

# check that this is indeed a walrus topology file
check = topfile.readline().rstrip()
print "file type = " + check
if check != 'walrus topology file':
        sys.exit("Error: walrus topology file is expected for argument 1")

# get the number of elements
check = topfile.readline().rstrip().split()
num_elem = int(check[1])
print "num_elements = " + str(num_elem)

topfile.readline() # skip "num_surface_elements" line
topfile.readline() # skip "num_interior_elements" line
topfile.readline() # skip "surface elements:" line

check = topfile.readline().rstrip()
if check != 'interior elements:':
        sys.exit("Error: Line 5 of topology file should read 'interior elements:'. Topology file is either badly constructed, or has already been run through later scripts...")

# read in all the elements
elements = topfile.readlines()

topfile.close()


# create an empty list to which the faces can be added
face_list = []

# Get the node numbers for each element and use these
# to construct the element's 16 faces. Order each face so
# that the node indices appear in ascending order (to
# make it easier to quickly check if two faces are the same)
print "Counting occurrences of each face..."
i = 0
for el in elements:
	a = el.split()
	for j in range(10):
		a[j] = int(a[j])

	# linear tetrahedron has 4 faces
#	add_face(face_list, a[0], a[1], a[2], i)
#	add_face(face_list, a[0], a[1], a[3], i)
#	add_face(face_list, a[0], a[2], a[3], i)
#	add_face(face_list, a[1], a[2], a[3], i)

	# quadratic tetrahedron has 16 faces
	add_face(face_list, a[3], a[6], a[8], i)
	add_face(face_list, a[0], a[4], a[6], i)
	add_face(face_list, a[4], a[6], a[8], i)
	add_face(face_list, a[1], a[4], a[8], i)

	add_face(face_list, a[0], a[4], a[5], i)
	add_face(face_list, a[1], a[4], a[7], i)
	add_face(face_list, a[4], a[5], a[7], i)
	add_face(face_list, a[2], a[5], a[7], i)

	add_face(face_list, a[3], a[6], a[9], i)
	add_face(face_list, a[0], a[5], a[6], i)
	add_face(face_list, a[5], a[6], a[9], i)
	add_face(face_list, a[2], a[5], a[9], i)

	add_face(face_list, a[3], a[8], a[9], i)
	add_face(face_list, a[1], a[7], a[8], i)
	add_face(face_list, a[7], a[8], a[9], i)
	add_face(face_list, a[2], a[7], a[9], i)

	i += 1

print "done."

print "Culling duplicates..."
for i in range(len(face_list)-1, -1, -1):
	face = face_list[i]
	if face[3] != 1:
		face_list.remove(face)
print "done."

print "Writing surface faces to file..."
surffile = open(sys.argv[2], "w")
surffile.write('walrus surface file\n')
surffile.write('num_surface_faces ' + str(len(face_list)) + '\n')
surffile.write('faces:\n')
for face in face_list:
	surffile.write(str(face[4]) + ' ' + str(face[0]) + ' ' + str(face[1]) + ' ' + str(face[2]) + '\n')
surffile.close()
print "Done. --> extract_surface_from_topology.py"
