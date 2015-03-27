#!/usr/bin/env python

import os, sys
import math

if len(sys.argv) != 3:
	sys.exit("Usage: python pdb_remove_type [INPUT .PDB FILE] [RADIUS THRESHOLD]")

ext = sys.argv[1][-4:]
base = sys.argv[1][0:-4]

if ext != ".pdb":
	sys.exit("Input file must be a '.pdb' file")

radius = float(sys.argv[2])
inputpdb = open(sys.argv[1], "r")
out_file = open(base + "_lessvolume" + ext, "w")
lines = inputpdb.readlines()
k = 0
for i in lines:
	if i[0:3] == "TER" and k == 1:
		i = i[0:7] + str(int(j[7:11]) + 1) + i[11:]
		out_file.write(i)  
		k = 0               
	elif i[0:4] != "ATOM" and i[0:6] != "HETATM": 	
		out_file.write(i)
	elif math.sqrt(math.pow(float(i[30:38].strip()), 2) + math.pow(float(i[38:46].strip()), 2)) < radius:	
		out_file.write(i)
		j = i
		k = 1
