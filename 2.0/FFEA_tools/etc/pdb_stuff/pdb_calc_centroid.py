#!/usr/bin/env python

import os, sys
import math

if len(sys.argv) != 2:
	sys.exit("Usage: python pdb_remove_type [INPUT .PDB FILE]")

inputpdb = open(sys.argv[1], "r")

if sys.argv[1][-4:] != ".pdb":
	sys.exit("Input file must be a '.pdb' file")

lines = inputpdb.readlines()
centroid = [0.0,0.0,0.0]
num_atoms = 0
for i in lines:
	if i[0:4] == "ATOM":
		centroid[0] += float(i[30:38].strip())
		centroid[1] += float(i[38:46].strip())
		centroid[2] += float(i[46:52].strip())
		num_atoms = num_atoms + 1

for i in range(3):
	centroid[i] = centroid[i] * 1.0/num_atoms

print "Num atoms = " + str(num_atoms) + "Centroid = ", centroid

