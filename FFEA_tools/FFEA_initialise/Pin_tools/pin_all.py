import sys, os

if len(sys.argv) != 3:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .node file] [OUTPUT .pin file]")

innodefname = sys.argv[1]
outpinfname = sys.argv[2]

# Get num_nodes
fin = open(innodefname, "r")
fin.readline()
num_nodes = int(fin.readline().split()[1])
fin.close()

# Pin all nodes
fout = open(outpinfname, "w")
fout.write("ffea pinned nodes file\n")
fout.write("num_pinned_nodes %d\n" % (num_nodes))
fout.write("pinned nodes:\n")
for i in range(num_nodes):
	fout.write(str(i) + "\n")
fout.close()
