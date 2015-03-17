#!/usr/bin/env python

import os, sys
import math

if len(sys.argv) != 4:
	sys.exit("Usage: python pdb_remove_type [INPUT .PDB FILE] [MODEL LOWER THRESHOLD] [MODEL UPPER THRESHOLD]")

ext = sys.argv[1][-4:]
base = sys.argv[1][0:-4]

if ext != ".pdb":
	sys.exit("Input file must be a '.pdb' file")

inputpdb = open(sys.argv[1], "r")
out_file = open(base + "_lessmodels" + ext, "w")
lines = inputpdb.readlines()
check = 2
for i in lines:
	if check == 2:
		if i[0:5] != "MODEL":
			out_file.write(i)
		else:
			check = 0

	if check == 0:
		if i[0:5] == "MODEL" and int(i[12:14]) == int(sys.argv[2]):
			check = 1
			out_file.write(i)
			continue
	else:
		if i[0:5] == "MODEL" and int(i[12:14]) == int(sys.argv[3]):
			check = 0	
		else:
			out_file.write(i)
			
