import sys, os
import math
import random

def read_pdb_frame(file):
        while True:
                line = file.readline()
                if line == "":
                        break
                line = line.replace("-", " -")
                if "REMARK" in line:
                        continue
                if "MODEL" in line:
                        continue
                if "ATOM" in line:
                        frame = []
                        while True:
                                if "ENDMDL" in line:
                                        break
                                if "TER" in line:
                                        break
                                if "END" in line:
                                        break
				# ignore the first 31 characters of the line
				line = line[31:]
				sline = line.split()

                                frame.append([float(sline[0]), float(sline[1]), float(sline[2])])
                                line = (file.readline())
                                line = line.replace("-", " -")
                        return len(frame), frame
        return -1, None

if len(sys.argv) != 4:
	sys.exit("Usage: python extract_nth_pdb_frame.py [INPUT PDB] [FRAME NUMBER TO OUTPUT] [OUTPUT PDB]")

pdb_in = open(sys.argv[1], "r")
N = int(sys.argv[2])

for i in xrange(N):
	num_atoms, frame = read_pdb_frame(pdb_in)

pdb_in.close()
pdb_out = open(sys.argv[3], "w")
for pos in frame:
	ridiculous_pos_x = ("%.3f" % (pos[0])).rjust(12, " ")
	ridiculous_pos_y = ("%.3f" % (pos[1])).rjust(8, " ")
	ridiculous_pos_z = ("%.3f" % (pos[2])).rjust(8, " ")

	stupid_number_format = str(i).rjust(7, " ")
	pdb_out.write("ATOM" + stupid_number_format + "  C   GLY     1" + ridiculous_pos_x + ridiculous_pos_y + ridiculous_pos_z + "\n")
	i += 1

pdb_out.write("TER\nENDMDL\n")

pdb_out.close()
