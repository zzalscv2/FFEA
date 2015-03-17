import sys

if len(sys.argv) != 5:
	sys.exit("Usage: python tetgen_to_walrus.py [INPUT .NODE FILE] [INPUT .ELE FILE] [OUTPUT NODE FILE] [OUTPUT TOPOLOGY FILE]")

inputnode = open(sys.argv[1],"r")
inputelem = open(sys.argv[2],"r")
outputnode = open(sys.argv[3], "w")
outputelem = open(sys.argv[4], "w")

# Get number of nodes
line = inputnode.readline().split()
num_nodes = int(line[0])
print "num_nodes = ", num_nodes

# Create node file
outputnode.write('walrus node file' + '\n')
outputnode.write('num_nodes ' + str(num_nodes) + '\n')
outputnode.write('num_surface_nodes ?\n')
outputnode.write('num_interior_nodes ?\n')
outputnode.write('surface nodes:\n')
outputnode.write('interior nodes:\n')

for i in range(num_nodes):
	line = inputnode.readline().split()
	outputnode.write(line[1] + ' ' + line[2] + ' ' + line[3] + '\n')

inputnode.close()
outputnode.close()

# Get number of elements
line = inputelem.readline().split()
num_elem = int(line[0])
print "num_elem = ", num_elem

# Create element file
outputelem.write('walrus topology file' + '\n')
outputelem.write('num_elements ' + str(num_elem) + '\n')
outputelem.write('num_surface_elements ?\n')
outputelem.write('num_interior_elements ?\n')
outputelem.write('surface elements:\n')
outputelem.write('interior elements:\n')
for i in range(num_elem):
	line = inputelem.readline().split()
	outputelem.write(line[1] + ' ' + line[2] + ' ' + line[3] + ' ' + line[4] + '\n')

inputelem.close()
outputelem.close()

print "tetgen_to_walrus.py -> Done."
