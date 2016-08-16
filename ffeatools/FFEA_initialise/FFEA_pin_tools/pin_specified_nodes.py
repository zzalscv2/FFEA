import sys, os

if len(sys.argv) < 2:
	sys.exit("Usage: python " + sys.argv[0] + " [OUTPUT .pin fname] [List of nodes to pin]")

outpinfname = sys.argv[1]
nodes_to_pin = []

for i in range(2, len(sys.argv)):
	nodes_to_pin.append(sys.argv[i])

# Pin all nodes
fout = open(outpinfname, "w")
fout.write("ffea pinned nodes file\n")
fout.write("num_pinned_nodes %d\n" % (len(nodes_to_pin)))
fout.write("pinned nodes:\n")
for i in range(len(nodes_to_pin)):
	fout.write(nodes_to_pin[i] + "\n")
fout.close()	
