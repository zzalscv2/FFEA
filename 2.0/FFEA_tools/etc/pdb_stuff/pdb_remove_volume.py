#!/usr/bin/env python

import os, sys

if len(sys.argv) != 8:
	sys.exit("Usage: python pdb_remove_type [INPUT .PDB FILE] [X_MAX] [X_MIN] [Y_MAX] [Y_MIN] [Z_MAX] [Z_MIN]")

ext = sys.argv[1][-4:]
base = sys.argv[1][0:-4]

if ext != ".pdb":
	sys.exit("Input file must be a '.pdb' file")

x_max = float(sys.argv[2])
x_min = float(sys.argv[3])
y_max = float(sys.argv[4])
y_min = float(sys.argv[5])
z_max = float(sys.argv[6])
z_min = float(sys.argv[7])
inputpdb = open(sys.argv[1], "r")
out_file = open(base + "_shrinked" + ext, "w")
lines = inputpdb.readlines()
k = 0
for i in lines:
	if i[0:3] == "TER" and k == 1:
		i = i[0:7] + str(int(j[7:11]) + 1) + i[11:]
		out_file.write(i)  
		k = 0               
	elif i[0:4] != "ATOM" and i[0:6] != "HETATM": 	
		out_file.write(i)
	elif float(i[30:38].strip()) > x_min and float(i[30:38].strip()) < x_max and float(i[38:46].strip()) > y_min and float(i[38:46].strip()) < y_max and float(i[46:54].strip()) > z_min and float(i[46:54].strip()) < z_max:	
		out_file.write(i)
		j = i
		k = 1
