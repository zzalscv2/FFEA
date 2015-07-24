#!/usr/bin/env python

import os, sys
import math

if len(sys.argv) != 3:
	sys.exit("Usage: python pdb_clean_for_vmd [INPUT .PDB FILE] [OUTPUT .PDB NAME]")

extin = sys.argv[1][-4:]
basein = sys.argv[1][0:-4]
extout = sys.argv[2][-4:]
baseout = sys.argv[2][0:-4]

if extin != ".pdb":
	sys.exit("Input file must be a '.pdb' file")

if extout != ".pdb":
	baseout = sys.argv[2]
	extout = ".pdb"

inputpdb = open(sys.argv[1], "r")
out_file = open(baseout + extout, "w")
lines = inputpdb.readlines()
check = 0
for i in lines:
	if check == 0:
		if i[0:4] == "COMP" or i[0:4] == "MODE" or i[0:4] == "HETA":
			out_file.write(i)
		elif i[0:4] == "ATOM":
			j = i
			check = 1
			out_file.write(i)
	else:
		if i[0:4] == "ATOM":
			out_file.write(i)
			j = i
		elif i[0:4] == "TER ":
			out_file.write(i)
			check = 0
		else:
			check = 0
			out_file.write("TER " + j[7:11] + "     " + j[17:26])
			out_file.write(i)

		


