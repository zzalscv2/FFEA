import sys, os

def PDB_get_num_nodes(fname):
	num_nodes = 0
	infile = open(fname, "r")
	while True:
		while "MODEL" not in infile.readline():
			continue
		while "TER" not in infile.readline():
			num_nodes += 1
		break

	return num_nodes

if len(sys.argv) != 4:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .pdb fname] [OUTPUT FFEA traj fname] [FFEA scale]")

infname = sys.argv[1]
outfname = sys.argv[2]
scale = float(sys.argv[3])

# Get number of nodes in pdb
num_nodes = PDB_get_num_nodes(infname)

# Open files
infile = open(infname, "r")
outfile = open(outfname, "w")

# Get number of nodes in 
# Outfile Header Stuff
outfile.write("FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs 1\nNumber of Conformations 1\nBlob 0: Conformation 0 Nodes %d\n\n" % (num_nodes))

step = -1

while True:
	line = infile.readline()
	if line == "":
		outfile.write("*")
		break

	while "MODEL" not in line:
		line = infile.readline()
		continue

	step += 1
	outfile.write("*\nBlob 0, Conformation 0, step %d\nDYNAMIC\n" % (step))
	line = infile.readline()
	while "TER" not in line:
		outfile.write("%6.3e %6.3e %6.3e " % (float(line[30:38]) * scale, float(line[38:46]) * scale, float(line[46:52]) * scale))	
		for i in range(7):
			outfile.write("%6.3e " % (0.0))
		outfile.write("\n")
		line = infile.readline()

	line = infile.readline()
	while "ENDMDL" not in line:
		line = infile.readline()
		continue
	
