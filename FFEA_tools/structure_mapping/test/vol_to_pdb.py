import sys, os
import math
import Vectors

if len(sys.argv) != 3:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .vol file] [OUTPUT .pdb file]")

infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")
outfile.write("COMPND   1 CHAIN: A\n")
while True:
	line = infile.readline()
	if line.strip() == "points":
		num_nodes = int(infile.readline())
		for i in range(num_nodes):
			sline = infile.readline().split()
			outfile.write("ATOM")
			outfile.write(str(i + 1).rjust(7))
			outfile.write("  C")
			outfile.write(str("ILE").rjust(6))
			outfile.write(str("A1").rjust(6))
			outfile.write("    ")
			for j in range(3):
				outfile.write("%8.3f" % float(sline[j]))
			outfile.write("\n")
		break

infile.close()
outfile.close()		
	
