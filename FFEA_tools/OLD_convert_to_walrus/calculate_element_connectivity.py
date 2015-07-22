import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python calculate_element_connectivity.py [INPUT TOPOLOGY FILE] [OUTPUT ELEMENT CONNECTIVITY FILE (VDW)]")

print "Running: calculate_element_connectivity.py"

top = open(sys.argv[1], "r")
vdw = open(sys.argv[2], "w")


# check that we have a walrus topology file
check = top.readline().rstrip()
print "file type = " + check
if check != 'walrus topology file':
        sys.exit("Error: walrus topology file is expected for argument 2")

# get total number of elements
check = top.readline().rstrip().split()
num_elements = int(check[1])

# get number of surface elements
check = top.readline().rstrip().split()
num_surface_elements = int(check[1])

# get number of interior elements
check = top.readline().rstrip().split()
num_interior_elements = int(check[1])

# check for 'surface elements:' line
check = top.readline().rstrip()
if check != 'surface elements:':
        sys.exit("Error: Line 5 of topology file should read 'surface elements:")

# put all the element node indices into a big list
elements = []
for i in range(num_surface_elements):
        elements.append(top.readline())

# check for 'interior elements:' line
check = top.readline().rstrip()
if check != 'interior elements:':
        sys.exit("Error: expecting 'interior elements:' line")

top.close()

print "Calculating which elements share a node..."
element_connectivity = [[] for i in range (num_surface_elements)]

for i in range(num_surface_elements):
	nodes_i = elements[i].split()
	for j in range(i, num_surface_elements):
		nodes_j = elements[j].split()
		yes = 0
		for ni in range(10):
			for nj in range(ni,10):
				if int(nodes_i[ni]) == int(nodes_j[nj]):
					yes = 1
					break
			if yes == 1:
				element_connectivity[i].append(j)
				break


vdw.write('walrus surface element connectivity file\n')
for el in element_connectivity:
	el.sort()
	num = len(el)
	el = ' '.join(map(str, el))
	vdw.write(str(num) + ": " + el + "\n")

vdw.close()

print "Done. --> calculate_element_connectivity.py"
