import sys

if len(sys.argv) != 2:
	sys.exit("Usage: python get_dimensions.py [INPUT NODE FILE]")

print "Running: get_dimensions.py"

infile = open(sys.argv[1], "r")

# check that infile is a node file
check = infile.readline().rstrip()
print "file type = " + check
if check != 'walrus node file':
        sys.exit("Error: walrus node file is expected for argument 1")

# get total number of nodes
check = infile.readline().rstrip().split()
num_nodes = int(check[1])

# get number of surface nodes
check = infile.readline().rstrip().split()
num_surface_nodes = int(check[1])

# get number of interior nodes
check = infile.readline().rstrip().split()
num_interior_nodes = int(check[1])

# check for 'surface nodes:' line
check = infile.readline().rstrip()
if check != 'surface nodes:':
        sys.exit("Error: Line 5 of node file should read 'surface nodes:")

# put all the node positions into a big list

x_min = float('inf')
x_max = -float('inf')
y_min = float('inf')
y_max = -float('inf')
z_min = float('inf')
z_max = -float('inf')
for i in range(num_surface_nodes):
	a = infile.readline().split()
	x = float(a[0])
	y = float(a[1])
	z = float(a[2])
	if x < x_min:
		x_min = x
	elif x > x_max:
		x_max = x
	if y < y_min:
		y_min = y
	elif y > y_max:
		y_max = y
	if z < z_min:
		z_min = z
	elif z > z_max:
		z_max = z

# check for 'interior nodes:' line
check = infile.readline().rstrip()
if check != 'interior nodes:':
        sys.exit("Error: expecting 'interior nodes:")

# put the remaining node positions into the big list
for i in range(num_interior_nodes):
	a = infile.readline().split()
	x = float(a[0])
	y = float(a[1])
	z = float(a[2])
	if x < x_min:
		x_min = x
	elif x > x_max:
		x_max = x
	if y < y_min:
		y_min = y
	elif y > y_max:
		y_max = y
	if z < z_min:
		z_min = z
	elif z > z_max:
		z_max = z


print "x_min = " + str(x_min)
print "x_max = " + str(x_max)
print "y_min = " + str(y_min)
print "y_max = " + str(y_max)
print "z_min = " + str(z_min)
print "z_max = " + str(z_max)

print "size_x = " + str(x_max - x_min)
print "size_y = " + str(y_max - y_min)
print "size_z = " + str(z_max - z_min)


infile.close()

print "Done. --> get_dimensions.py"
