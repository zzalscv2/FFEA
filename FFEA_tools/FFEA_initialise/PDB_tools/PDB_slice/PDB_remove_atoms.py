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

printing = 1
num_lost_atoms = 0
residue_diff = 0
for i in range(len(lines)):

	# Dead line
	if lines[i] == "" or lines[i] == "\n":
		continue
	sline = lines[i].split()

	if sline[0].strip() == "ATOM":

		# Do we need to stop/start printing
		if int(sline[1].strip()) == first_atom:
			printing = 0
			start_residue = int(sline[4][1:])
		elif int(sline[1].strip()) == last_atom + 1:
			printing = 1

		# What to do in both cases
		if printing == 0:
			num_lost_atoms += 1
		else:
			lines[i] = lines[i].replace(sline[1], str(int(sline[1]) - num_lost_atoms).rjust(len(sline[1])))
			outfile.write(lines[i])

	elif sline[0].strip() == "TER":
		outfile.write(lines[i])
	else:
		continue

#sline = lines[end].split()
#outfile.write("TER   " + str(int(sline[1].strip()) + 1) + "      " + sline[3].strip() + " " + sline[4].strip() + "\n")
outfile.close()
