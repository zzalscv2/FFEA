import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python make_default_vdw.py [INPUT SURFACE FILE] [OUTPUT VDW FILE]")

with open(sys.argv[1], "r") as surf:
	surf.readline()
	line = surf.readline().split()
	num_faces = int(line[1])

vdw = open(sys.argv[2], "w")
vdw.write("ffea vdw file\n")
vdw.write("num_faces " + str(num_faces) + "\n")
vdw.write("vdw params:\n")
for i in range(num_faces):
	vdw.write("-1\n")
vdw.close()
