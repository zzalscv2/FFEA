import sys

if len(sys.argv) != 4:
	sys.exit("Usage: python DAT_to_walrus.py [VOL FILE] [OUTPUT NODE FILE] [OUTPUT TOPOLOGY FILE]")

print "Running: DAT_to_walrus.py"

inputfile = open(sys.argv[1],"r")
nodefile = open(sys.argv[2], "w")
topfile = open(sys.argv[3], "w")

num_nodes, num_elem = inputfile.readline().split()
num_nodes = int(num_nodes)
num_elem = int(num_elem)
print "num_nodes =", num_nodes
print "num_elements =", num_elem

# write nodes file
nodefile.write('walrus node file\n')
nodefile.write('num_nodes ' + str(num_nodes) + '\n')
nodefile.write('num_surface_nodes ?\n')
nodefile.write('num_interior_nodes ?\n')

print "Getting nodes..."
nodefile.write('surface nodes:\n')
nodefile.write('interior nodes:\n')
for i in range(num_nodes):
	line = inputfile.readline()
        a = line.split()
	nodefile.write(a[1] + ' ' + a[2] + ' ' + a[3] + '\n')
nodefile.close()


# write header for topology file
topfile.write('walrus topology file\n')
topfile.write('num_elements ' + str(num_elem) + '\n')
topfile.write('num_surface_elements ?\n')
topfile.write('num_interior_elements ?\n')

# get the node numbers for each element
print "Getting elements..."
topfile.write('surface elements:\n')
topfile.write('interior elements:\n')
for i in range(num_elem):
	line = inputfile.readline()
	a = line.split()
	topfile.write(str(int(a[1]) - 1) + ' ' + str(int(a[2]) - 1) + ' ' + str(int(a[3]) - 1) + ' ' + str(int(a[4]) - 1) + '\n')
topfile.close()
print "done"

inputfile.close()
print "Done. --> netgen_to_walrus.py"
