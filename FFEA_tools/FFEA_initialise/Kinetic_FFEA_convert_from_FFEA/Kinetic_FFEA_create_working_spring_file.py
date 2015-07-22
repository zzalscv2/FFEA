import sys, os

if len(sys.argv) < 5:
	sys.exit("Usage: python " + sys.argv[0] + " [spring_file_fname] [spring_constant] [tolerance_length] [List of pairs of nodes]")

spring_file_fname = sys.argv[1]
k = float(sys.argv[2])
l = float(sys.argv[3])
node_pairs = []
for i in range(4, len(sys.argv), 2):
	node_pairs.append([sys.argv[i], sys.argv[i + 1]])

fout = open(spring_file_fname, "w")
fout.write("ffea spring file\n")
fout.write("num_springs %d\n" % (len(node_pairs)))
fout.write("springs:\n")
for node_pair in node_pairs:
	fout.write("0 0 %d 1 0 %d %e %e\n" % (int(node_pair[0]), int(node_pair[1]), k, l))

fout.close()
