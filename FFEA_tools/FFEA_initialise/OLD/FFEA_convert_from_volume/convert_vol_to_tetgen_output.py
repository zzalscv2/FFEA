import os, sys

if len(sys.argv) != 4:
	sys.exit("Usage: python convert_vol_to_tetgen_output.py [INPUT VOL FILE] [OUTPUT .NODE FILE] [OUTPUT .ELE FILE]")

inputvol = open(sys.argv[1], "r")
outputnode = open(sys.argv[2],"w")
outputelem = open(sys.argv[3],"w")

lines = inputvol.readlines()
i = lines.index("volumeelements\n")
num_elem = int(lines[i+1])
elem = []
for j in range(i+2, i + 2 + num_elem, 1):
	sline = [int(ei) - 1 for ei in lines[j].split()]
	s = str(sline[2]) + " " + str(sline[3]) + " " + str(sline[4]) + " " + str(sline[5])
	elem.append(s)

i = lines.index("points\n")
num_nodes = int(lines[i+1])
nodes = []
for j in range(i+2, i + 2 + num_nodes, 1):
	sline = lines[j].split()
	s = str(sline[0]) + " " + str(sline[1]) + " " + str(sline[2])
	nodes.append(s)


outputnode.write(str(num_nodes) + " 3 0 0\n")
for n in range(len(nodes)):
	outputnode.write(str(n) + " " + nodes[n] + "\n")

outputelem.write(str(num_elem) + " 4 0\n")
for e in range(len(elem)):
	outputelem.write(str(e) + " " + elem[e] + "\n")

print "convert_vol_to_tetgen_output.py -> Done."
