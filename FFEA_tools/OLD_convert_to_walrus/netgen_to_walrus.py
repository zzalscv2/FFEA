import sys

if len(sys.argv) != 4:
	sys.exit("Usage: python netgen_to_walrus.py [VOL FILE] [OUTPUT NODE FILE] [OUTPUT TOPOLOGY FILE]")

print "Running: netgen_to_walrus.py"

inputfile = open(sys.argv[1],"r")
nodefile = open(sys.argv[2], "w")
topfile = open(sys.argv[3], "w")

# read through file until line that says "volumeelements"
while inputfile.readline().find('volumeelements') == -1:
		pass

# get number of elements
num_elem = int(inputfile.readline())
topfile.write('walrus topology file\n')
topfile.write('num_elements ' + str(num_elem) + '\n')
topfile.write('num_surface_elements ?\n')
topfile.write('num_interior_elements ?\n')
print str(num_elem) + ' elements'

# get the node numbers for each element
print "Getting elements..."
topfile.write('surface elements:\n')
topfile.write('interior elements:\n')
for i in range(num_elem):
	line = inputfile.readline()
	a = line.split()
	topfile.write(str(int(a[2]) - 1) + ' ' + str(int(a[3]) - 1) + ' ' + str(int(a[4]) - 1) + ' ' + str(int(a[5]) - 1) + '\n')
topfile.close()
print "done"

# read through file until line that says "points"
while inputfile.readline().find('points') == -1:
		pass

# get number of elements
num_nodes = int(inputfile.readline())
nodefile.write('walrus node file\n')
nodefile.write('num_nodes ' + str(num_nodes) + '\n')
nodefile.write('num_surface_nodes ?\n')
nodefile.write('num_interior_nodes ?\n')
print str(num_nodes) + ' nodes'

print "Getting nodes..."
nodefile.write('surface nodes:\n')
nodefile.write('interior nodes:\n')
for i in range(num_nodes):
	line = inputfile.readline()
        a = line.split()
	nodefile.write(a[0] + ' ' + a[1] + ' ' + a[2] + '\n')
nodefile.close()

inputfile.close()
print "Done. --> netgen_to_walrus.py"
