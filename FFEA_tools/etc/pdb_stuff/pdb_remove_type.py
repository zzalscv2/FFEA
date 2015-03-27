#!/usr/bin/env python

import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python pdb_remove_type [INPUT .PDB FILE] [PDB TYPE]")

ext = sys.argv[1][-4:]
base = sys.argv[1][0:-4]

if ext != ".pdb":
	sys.exit("Input file must be a '.pdb' file")
inputpdb = open(sys.argv[1], "r")
out_file = open(base + "_no" + sys.argv[2] + "_" + ext, "w")
lines = inputpdb.readlines()

for i in lines:
	if i[0:len(sys.argv[2])] != sys.argv[2]:	
		out_file.write(i)
		
