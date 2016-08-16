import sys, os

if len(sys.argv) != 5:
	sys.exit("Usage python " + sys.argv[0] + " [INPUT .pdb fname] [OUTPUT .pdb fname] [from atom] [to atom]")

# Get args
infname = sys.argv[1]
outfname = sys.argv[2]
first_atom = int(sys.argv[3])
last_atom = int(sys.argv[4])

# Open both file for immediate appending
infile = open(infname, "r")
outfile = open(outfname, "w")
outfile.write("REMARK   1 Created by '" + os.path.basename(sys.argv[0]) + "'\n")

# Run through all lines of infile, printing only required atoms, ignoring rest
lines = infile.readlines()
infile.close()

for i in range(len(lines)):
	if lines[i] == "" or lines[i] == "\n":
		continue
	sline = lines[i].split()
	if sline[0].strip() != "ATOM":
		continue

	if int(sline[1].strip()) == first_atom:
		start_line = i
		break

for i in range(start_line, len(lines)):
	if lines[i] == "" or lines[i] == "\n":
		continue

	sline = lines[i].split()
	if sline[0] == "ATOM" and int(sline[1]) <= last_atom:
		outfile.write(lines[i])

	if int(sline[1]) == last_atom:
		end_line = i
		break

sline = lines[end_line].split()
outfile.write("TER   " + str(int(sline[1].strip()) + 1) + "      " + sline[3].strip() + " " + sline[4].strip() + "\n")
outfile.close()
