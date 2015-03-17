import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python convert_pdb_to_walrus.py [INPUT PDB FILE] [OUTPUT WALRUS TRAJ FILE]")

# Get the number of nodes in the pdb (extremely inefficiently, but who gives a shit? pdb files have a mind bendingly stupid format.)
pdb = open(sys.argv[1], "r")
num_nodes = -1
while True:
	line = pdb.readline().split()
	if "REMARK" in line[0] or "MODEL" in line[0]:
		if num_nodes == -1:
			continue
		else:
			break
	if "ATOM" in line[0]:
		num_nodes = int(line[1])
			
pdb.close()
num_nodes += 1
print "num_nodes =", num_nodes

pdb = open(sys.argv[1], "r")
walrus_traj = open(sys.argv[2], "w")
frame = 0
while True:
	line = pdb.readline()
	if line == "":
		break
	line = line.replace("-", " -")
	line = line.split()
	if "REMARK" in line[0]:
		continue
	if "MODEL" in line[0]:
		continue
	if "ATOM" in line[0]:
		print "Converting frame ", frame
		walrus_traj.write("*\nBlob 0, step " + str(frame) + "\nDYNAMIC\n" + str(num_nodes) + "\n")
		while "ENDMDL" not in line[0]:
			walrus_traj.write(str(float(line[5]) * 1e-9) + " " + str(float(line[6]) * 1e-9) + " " + str(float(line[7]) * 1e-9) + " 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00\n")
			line = (pdb.readline())
			line = line.replace("-", " -")
			line = line.split()
		frame += 1

pdb.close()
walrus_traj.close()
