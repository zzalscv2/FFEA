#!/usr/bin/env python

import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python pdb_chain_convert [ORIGINAL .PDB FILE] [SECTION TO CHANGE .PDB FILE]")

ext1 = sys.argv[1][-4:]
base1 = sys.argv[1][0:-4]
ext2 = sys.argv[2][-4:]
base2 = sys.argv[2][0:-4]
if ext1 != ".pdb" or ext2 != ".pdb":
	sys.exit("Both input files must be '.pdb' files")

inputpdb = open(sys.argv[1], "r")
changepdb = open(sys.argv[2], "r")
out_file = open(base1 + "_chain_converted" + ext1, "w")
lines1 = inputpdb.readlines()
lines2 = changepdb.readlines()
num_chains = 0
chain = []
for i in lines1:
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

convert_to = "D"
for i in lines1:
	if i[0:4] == "ATOM":		
		for j in lines2: 
			if j[0:4] == "ATOM":
				if i[30:38] == j[30:38] and i[38:46] == j[38:46] and i[46:54] == j[46:54]:
					k = i[0:21] + convert_to + i[22:]
					out_file.write(k)
					break
	elif i[0:6] == "HETATM":
		for j in lines2: 
			if j[0:4] == "HETATM":
				if i[30:38] == j[30:38] and i[38:46] == j[38:46] and i[46:54] == j[46:54]:
					k = i[0:21] + convert_to + i[22:]
					out_file.write(k)
					break
	else:
		out_file.write(i)








