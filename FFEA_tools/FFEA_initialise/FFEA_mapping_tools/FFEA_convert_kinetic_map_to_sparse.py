import sys, os

if len(sys.argv) != 3:
	sys.exit("Usage: python " + os.path.basename(sys.argv[0]) + " [INPUT FFEA kinetic map fname] [OUTPUT FFEA kinetic map fname]")

# Get args
map_fname = sys.argv[1]
outmap_fname = sys.argv[2]

# Open and read header info
inmap = open(map_fname, "r")
line = inmap.readline().strip()
if line != "FFEA Kinetic Conformation Mapping File (Dense)":
	sys.exit("Expected 'FFEA Kinetic Conformation Mapping File (Dense)' and got " + line + ". May not be the correct file type\n")

num_columns = int(inmap.readline().split()[1])
num_rows = int(inmap.readline().split()[1])
num_entries_start = int(inmap.readline().split()[1])

# Last line
inmap.readline()

# Define sparse containers
A = []
IA = []
JA = []
IA.append(0)

for i in range(num_rows):
	if i >= 100:
		if i % int(num_rows / 100.0) == 0:
			sys.stdout.write("\r%d%% read" % (int(i * 100.0 / num_rows)))
			sys.stdout.flush()
	sline = inmap.readline().split()
	for j in range(len(sline)):
		if float(sline[j]) != 0.0:
			A.append(sline[j])
			JA.append(j)

	IA.append(len(A))

inmap.close()

if (num_entries_start != len(A)):
	sys.exit("Error. Specified num_entries not equal to num_read")
sys.stdout.write("\r100% read\n")
sys.stdout.flush()

outmap = open(outmap_fname, "w")
outmap.write("FFEA Kinetic Conformation Mapping File (Sparse)\n")
outmap.write("num_nodes_from %d\n" % (num_columns))
outmap.write("num_nodes_to %d\n" % (num_rows))
outmap.write("num_entries %d\n" % (len(A)))
outmap.write("map:\n")
outmap.write("entries - ")
for i in A:
	outmap.write(i + " ")
outmap.write("\n")
outmap.write("key - ")
for i in IA:
	outmap.write(str(i) + " ")
outmap.write("\n")
outmap.write("columns - ")
for i in JA:
	outmap.write(str(i) + " ")
outmap.write("\n")
outmap.close()
outmap = open(outmap_fname, "r")
for i in range(5):
	outmap.readline()
