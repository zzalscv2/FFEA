#!/usr/bin/env python

import os, sys
import math

if len(sys.argv) != 4:
	sys.exit("Usage: python pdb_remove_type [INPUT .PDB FILE] [START ATOM TO REMOVE] [END ATOM TO REMOVE]")

ext = sys.argv[1][-4:]
base = sys.argv[1][0:-4]

if ext != ".pdb":
	sys.exit("Input file must be a '.pdb' file")

inputpdb = open(sys.argv[1], "r")
outfile = open(base + "_lessatoms" + ext, "w")
lines = inputpdb.readlines()

lower = int(sys.argv[2])
upper = int(sys.argv[3])
offset = upper - lower

check = 0;

for i in lines:
	if check == 0:
		if int(i[6:11]) == lower:
			check = 1
		else:
			outfile.write(i)
	elif check == 1:
		if int(i[6:11]) == upper:
			check = 2
	elif check == 2:
		new_i = i[0:4] + str(int(i[6:11]) - offset).rjust(7) + i[11:]
		outfile.write(new_i)
