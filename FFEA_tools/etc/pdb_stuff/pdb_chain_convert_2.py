#!/usr/bin/env python

import os, sys
import math

if len(sys.argv) != 3:
	sys.exit("Usage: python pdb_remove_type [INPUT .PDB FILE] [RADIUS]")

ext = sys.argv[1][-4:]
base = sys.argv[1][0:-4]

if ext != ".pdb":
	sys.exit("Input file must be a '.pdb' file")

inputpdb = open(sys.argv[1], "r")
out_file = open(base + "_converted" + ext, "w")
lines = inputpdb.readlines()
convert_to = "D"
for i in lines:
	if i[0:4] != "ATOM" and i[0:6] != "HETATM":
		out_file.write(i) 	
	else:
		radius = math.sqrt(math.pow(float(i[30:38].strip()), 2) + math.pow(float(i[38:46].strip()), 2))
		if radius > float(sys.argv[2]):	
			k = i[0:21] + convert_to + i[22:]
			out_file.write(k)
		else:
			out_file.write(i)
