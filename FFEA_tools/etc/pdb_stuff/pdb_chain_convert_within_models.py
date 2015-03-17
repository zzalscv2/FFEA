#!/usr/bin/env python

import os, sys

if len(sys.argv) != 4:
	sys.exit("Usage: python pdb_chain_convert_within_models [ORIGINAL .PDB FILE] [MODEL LOWER THRESHHOLD] [MODEL UPPER THRESHOLD]")

ext = sys.argv[1][-4:]
base = sys.argv[1][0:-4]
if ext != ".pdb":
	sys.exit("Input file must be a '.pdb' file")

inputpdb = open(sys.argv[1], "r")
out_file = open(base + "_chain_converted" + ext, "w")
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

print "Found " + str(num_chains) + " chains."
for i in chain:
	print "Chain " + i

check = 0 
while check == 0:
	convert_from = raw_input("Convert from?:")
	for i in chain:
		if convert_from == i:
			check = 1
	if check == 0:
		print "Please select a valid chain identifier."

check = 0 
while check == 0:
	convert_to = raw_input("Convert to?:")
	for i in chain:
		if convert_to == i:
			check = 1
	if check == 0:
		print "Please select a valid chain identifier."
			
check = 0
for i in lines:
	if check == 0:
		out_file.write(i)
		if i[0:5] == "MODEL" and int(i[12:14]) == int(sys.argv[2]):
			check = 1
			continue
	else:
		if i[0:4] == "ATOM" or i[0:6] == "HETATM":
			if i[21] == convert_from:
				out_file.write(i[0:21] + convert_to + i[22:])
			else:
				out_file.write(i)

		elif i[0:5] == "MODEL" and int(i[12:14]) == int(sys.argv[3]):
			out_file.write(i)
			check = 0	
		else:
			out_file.write(i)
	







