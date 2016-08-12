import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python convert_obj_to_off.py [INPUT OBJ FNAME] [OUTPUT OFF FNAME]")

obj = open(sys.argv[1], "r")
v = []
f = []
for line in obj.readlines():
	if line.split()[0] == "v":
		v.append(line.replace("v ",""))
	if line.split()[0] == "f":
		if "//" in line:
			line = line.replace("f ","").replace("//", " ")
			line = line.split()
			f.append(str(int(line[0])-1) + " " + str(int(line[2])-1) + " " + str(int(line[4])-1) + "\n")
		else:
			line = line.replace("f ", "")
			line = line.split()
			f.append(str(int(line[0])-1) + " " + str(int(line[1])-1) + " " + str(int(line[2])-1) + "\n")
obj.close()

off = open(sys.argv[2], "w")
off.write("OFF\n")
off.write(str(len(v)) + " " + str(len(f)) + " 0\n")
for vert in v:
	off.write(vert)
for face in f:
	off.write("3 " + face)
off.close()
