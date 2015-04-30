import sys

if len(sys.argv) != 4:
	sys.exit("Usage: python oliver_to_walrus.py [OLIVER FILE] [OUTPUT NODE FILE] [OUTPUT TOPOLOGY FILE]")

inputfile = open(sys.argv[1],"r")
nodefile = open(sys.argv[2], "w")
topfile = open(sys.argv[3], "w")

# read the first line and get node and element numbers
line = inputfile.readline()
a = line.split()
num_nodes = int(a[0])
num_elem = int(a[1])

print str(num_nodes) + ' nodes, ' + str(num_elem) + ' elements'
nodefile.write(str(num_nodes) + '\n')
topfile.write(str(num_elem) + '\n')
print "Converting Oliver file format..."

print "Getting nodes..."
for i in range(num_nodes):
	line = inputfile.readline()
        a = line.split()
	x = float(a[1]) * 1.0e-10
	y = float(a[2]) * 1.0e-10
	z = float(a[3]) * 1.0e-10
	nodefile.write(str(x) + ' ' + str(y) + ' ' + str(z) + '\n')
nodefile.close()

# remember to subtract 1 from node indices as numbering in walrus starts from 0
print "Getting elements..."
for i in range(num_elem):
	line = inputfile.readline()
        a = line.split()
	topfile.write(str(int(a[1]) - 1) + ' ' + str(int(a[2]) - 1) + ' ' + str(int(a[3]) - 1) + ' ' + str(int(a[4]) - 1) + '\n')
topfile.close()

inputfile.close()

print "Done"
