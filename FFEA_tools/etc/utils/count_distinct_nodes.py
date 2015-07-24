import sys

def add_node(nl, pn):
	for n in nl:
		if pn == n[0]:
			return
	nl.append([pn]);

if len(sys.argv) != 2:
	sys.exit("Usage: python count_distinct_nodes.py [SURFACE FILE]")

surffile = open(sys.argv[1], "r")

# get number of elements listed in surface file
num_elem = int(surffile.readline())

# create an empty list to which the distinct nodes can be added
node_list = []

print "Constructing distinct node list..."
for i in range(num_elem):
	line = surffile.readline()
	a = line.split()
	a[1] = int(a[1])
	a[2] = int(a[2])
	a[3] = int(a[3])
	add_node(node_list, a[1])
	add_node(node_list, a[2])
	add_node(node_list, a[3])
surffile.close()
print "done."

print "Sorting node list..."
node_list.sort();
print "done."

for node in node_list:
	print node

print "done."

print "Num distinct nodes = " + str(len(node_list))
