import sys

if len(sys.argv) != 4:
	sys.exit("Usage: python scale.py [INPUT NODE FILE] [OUTPUT NODE FILE] [SCALING FACTOR]")

print "Running: scale.py"

infile = open(sys.argv[1], "r")
scaling = float(sys.argv[3])

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
node = []
for i in range(num_surface_nodes):
        a = infile.readline().split()
        node.append([float(a[0]), float(a[1]), float(a[2])])

# check for 'interior nodes:' line
check = infile.readline().rstrip()
if check != 'interior nodes:':
        sys.exit("Error: expecting 'interior nodes:")

# put the remaining node positions into the big list
for i in range(num_interior_nodes):
        a = infile.readline().split()
        node.append([float(a[0]), float(a[1]), float(a[2])])

infile.close()

# create the new, scaled node file
outfile = open(sys.argv[2], "w")
outfile.write("walrus node file\n")
outfile.write("num_nodes " + str(num_nodes) + "\n")
outfile.write("num_surface_nodes " + str(num_surface_nodes) + "\n")
outfile.write("num_interior_nodes " + str(num_interior_nodes) + "\n")
outfile.write("surface nodes:\n")

count = 0
for n in node:
	nx = n[0] * scaling
	ny = n[1] * scaling
	nz = n[2] * scaling
	outfile.write(str(nx) + " " + str(ny) + " " + str(nz) + "\n")
	count += 1
	if count == num_surface_nodes:
		outfile.write("interior nodes:\n")

outfile.close()

print "Done. --> scale.py"
