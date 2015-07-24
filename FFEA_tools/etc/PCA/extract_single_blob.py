import os, sys

if len(sys.argv) != 4:
	sys.exit("Usage: python extract_single_blob.py [INPUT FFEA TRAJECTORY FILE] [OUTPUT FFEA TRAJECTORY FILE] [BLOB NUMBER]")

infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")
blob_number = sys.argv[3]

infile.readline() #*
outfile.write("*\n")
check = 0

while True:
	line = infile.readline()
	line = line.strip().split()
	try:
		if line[0] == "*":
			outfile.write("*\n")
			continue
	except IndexError:
		print "Index Error signals eof. Done!"
		break

	if line[1][0] == blob_number:
		outfile.write(line[0] + " " + line[1] + " " + line[2] + " " + line[3] + "\n")
		check = check + 1
		if check == 100:
			check = 0
			print "Step " + line[3]
		line = infile.readline().strip()
		if line == "STATIC":
			outfile.write(line + "\n")
			continue
		else:
			outfile.write(line + "\n")
			num_nodes = int(infile.readline())
			outfile.write(str(num_nodes) + "\n")
			for i in range(num_nodes):
				outfile.write(infile.readline())
			continue
	if line[1][0] != blob_number:
		line = infile.readline().strip()
		if line == "STATIC":
			continue
		else:
			num_nodes = int(infile.readline())
			for i in range(num_nodes):
				infile.readline()
			continue

infile.close()
outfile.close()
