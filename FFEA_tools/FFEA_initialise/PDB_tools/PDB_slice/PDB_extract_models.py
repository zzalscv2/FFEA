import sys

if len(sys.argv) != 5:
	sys.exit("Usage: python PDB_extract_models.py [INPUT .pdb fname] [OUTPUT .pdb fname] [from model] [to model]")

# Get args
infname = sys.argv[1]
outfname = sys.argv[2]
from_model = int(sys.argv[3])
to_model = int(sys.argv[4])

# Open files and begin reading
fin = open(infname, "r")
fout = open(outfname, "w")

completed = 0
while completed == 0:
	line = fin.readline()
	if line[0:5] == "MODEL" and int(line.split()[1]) == from_model:
		fout.write(line)
		completed = 1

line = fin.readline()
while not (line[0:5] == "MODEL" and int(line.split()[1]) == to_model):
	fout.write(line)
	line = fin.readline()
