import sys, os
from Vectors import vector3

if len(sys.argv) != 6:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .vol file] [OUTPUT .vol file] [x shift] [y shift] [z shift]")


# Get args
infname = sys.argv[1]
outfname = sys.argv[2]
translation = vector3(float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))

fin = open(infname, "r")
fout = open(outfname, "w")

# Get to points bit
while True:
	line = fin.readline()
	fout.write(line)
	if line.strip() == "points":
		break

# Get num_points
line = fin.readline()
fout.write(line)
num_nodes = int(line.strip())

# Write translated nodes
for i in range(num_nodes):
	line = fin.readline()
	sline = line.split()
	node = vector3(float(sline[0]), float(sline[1]), float(sline[2]))
	node += translation
	fout.write("%6.3e %6.3e %6.3e\n" % (node.x, node.y, node.z))

fin.close()
fout.close()
