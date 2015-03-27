#!/usr/bin/env python

import os, sys

if len(sys.argv) != 2:
	sys.exit("Usage: python pdb_separate_chains [INPUT .PDB FILE]")

ext = sys.argv[1][-4:]
base = sys.argv[1][0:-4]

if ext != ".pdb":
	sys.exit("Input file must be a '.pdb' file")
inputpdb = open(sys.argv[1], "r")

lines = inputpdb.readlines()
num_chains = 0
chain = []
for i in lines:
	if i[0:6] == "COMPND":
		if i[11:16] == "CHAIN":
			count = 16
			while i[count] != ";":
				if i[count] == ":":
					count += 1
					continue
				elif i[count] == " ":
					count += 1
					continue
				elif i[count] == ",":
					count += 1
					continue
				else :
					chain.append(i[count])					
					num_chains += 1
					count += 1
			print "Num_chains: " + str(num_chains)
			break

if num_chains == 0:
	sys.exit("No chains explicitly specified. Please add line of the form 'COMPND *** CHAIN: A, B, C...;'")

out_file = []
for i in range(0, num_chains):
	out_file.append(open(base + "_chain_" + chain[i] + ext, "w"))

for i in lines:
	if i[0:6] == "COMPND":
		for j in range(0, num_chains):
			out_file[j].write("COMPND   1 CHAIN: " + chain[j] + ";\n")
	elif i[0:5] == "HELIX":  
		for j in range(0, num_chains):
			if i[19] == chain[j]:
				out_file[j].write(i)
	elif i[0:5] == "SHEET":
		 for j in range(0, num_chains):
			if i[21] == chain[j]:
				out_file[j].write(i)	
	elif i[0:6] == "MODEL":
		for j in range(0, num_chains):
			out_file[j].write(i)
	elif i[0:4] == "ATOM":
		 for j in range(0, num_chains):
			if i[21] == chain[j]:
				out_file[j].write(i)	
	elif i[0:3] == "TER":
		 for j in range(0, num_chains):
			if i[21] == chain[j]:
				out_file[j].write(i)
	elif i[0:6] == "HETATM":
		 for j in range(0, num_chains):
			if i[21] == chain[j]:
				out_file[j].write(i)
	if i[0:3] == "END":
		for j in range(0, num_chains):
			out_file[j].write(i)

		




