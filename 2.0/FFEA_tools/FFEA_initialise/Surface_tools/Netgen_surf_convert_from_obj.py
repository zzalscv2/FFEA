import os, sys

if len(sys.argv) != 3:
	sys.exit("Usage: python convert_obj_to_surf.py [INPUT OBJ FNAME] [OUTPUT SURF FNAME]")

obj = open(sys.argv[1], "r")
v = []
f = []
for line in obj.readlines():
	if line[0] != "v" and line[0] != "f":
		continue
	one_line = line.split()
	if one_line[0] == "v":
		v.append(line.replace("v ",""))
	if one_line[0] == "f":
		index = [0,0,0]
		index[0] = int(one_line[1].split("//")[0])
		index[1] = int(one_line[2].split("//")[0])
		index[2] = int(one_line[3].split("//")[0])
		for i in range(len(index)):
			if index[i] < 0:
				index[i] = -1 * index[i]

		face_str = str(index[0]) + " " + str(index[1]) + " " + str(index[2]) + "\n"
		f.append(face_str)
obj.close()

surf = open(sys.argv[2], "w")
surf.write("surfacemesh\n")
surf.write(str(len(v)) + "\n")
for vert in v:
	surf.write(vert)
surf.write(str(len(f)) + "\n")
for face in f:
	surf.write(face)
surf.close()
