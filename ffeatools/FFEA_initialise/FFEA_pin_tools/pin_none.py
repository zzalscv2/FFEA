import sys, os

if len(sys.argv) != 2:
	sys.exit("Usage: python " + sys.argv[0] + " [OUTPUT .pin file]")

outpinfname = sys.argv[1]

# Pin 0 nodes
fout = open(outpinfname, "w")
fout.write("ffea pinned nodes file\n")
fout.write("num_pinned_nodes 0\n")
fout.write("pinned nodes:\n")
fout.close()
