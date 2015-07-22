import sys

if len(sys.argv) != 5:
	sys.exit("Usage: python reset_vdw.py [INPUT VDW FILE] [OUTPUT VDW FILE] [FROM TYPE] [TO TYPE]")

inputvdw = open(sys.argv[1], "r")
outputvdw = open(sys.argv[2], "w")
fromtype = sys.argv[3]
totype = sys.argv[4]

#Initial crap
line = inputvdw.readline()
outputvdw.write(line)
line = inputvdw.readline()
num_faces = int(line.split()[1])
outputvdw.write(line)
line = inputvdw.readline()
outputvdw.write(line)

#Change vdw type
for i in range(num_faces):
	line = inputvdw.readline()
	if line.strip() == fromtype:
		outputvdw.write(totype + "\n")
	else:
		outputvdw.write(line)

inputvdw.close()
outputvdw.close()
			


