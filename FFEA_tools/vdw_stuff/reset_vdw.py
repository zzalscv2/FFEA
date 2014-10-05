import sys

if len(sys.argv) != 2:
	sys.exit("Usage: python reset_vdw.py [INPUT VDW FILE]")

inputvdw = open(sys.argv[1], "r")

line = inputvdw.readline()
num_faces = int(inputvdw.readline().split()[1])
inputvdw.close()
outputvdw = open(sys.argv[1], "w")
outputvdw.write("walrus vdw file\n")
outputvdw.write("num_faces " + str(num_faces) + "\n")
outputvdw.write("vdw params:\n")

for i in range(0, num_faces):
	outputvdw.write("-1\n")

