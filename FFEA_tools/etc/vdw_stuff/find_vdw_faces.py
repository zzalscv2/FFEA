import sys

if len(sys.argv) != 2:
	sys.exit("Usage: python find_vdw_faces.py [INPUT VDW FILE]")

inputvdw = open(sys.argv[1], "r")
lines = inputvdw.readlines()

i = 0
for line in lines:
	if line == "1\n" or line == "1":
		print i - 3
	i = i + 1
