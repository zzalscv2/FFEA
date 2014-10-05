import os, sys

if len(sys.argv) != 4:
	sys.exit("Usage: python convert_obj_to_FFEA_surf.py [INPUT OBJ FNAME] [OUTPUT FFEA SURF FNAME] [OUTPUT FFEA NODES FNAME]")

obj = open(sys.argv[1], "r")
v = []
f = []
for line in obj.readlines():
	if line == "\n":
		continue
	if line.split()[0] == "v":
		v.append(line.replace("v ",""))
	if line.split()[0] == "f":
		f.append(line.replace("f ",""))
obj.close()

surf = open(sys.argv[2], "w")
surf.write("walrus surface file\n")
surf.write("num_surface_faces " + str(len(f)) + "\n")
surf.write("faces:\n")
for face in f:
	face_split = face.split()
	index = [0,0,0]
	index[0] = face_split[0].split("//")[0]
	index[1] = face_split[1].split("//")[0]
	index[2] = face_split[2].split("//")[0]
	face_str = "0 " + str(int(index[0]) - 1) + " " + str(int(index[1]) - 1) + " " + str(int(index[2]) - 1) + "\n"
	surf.write(face_str)
surf.close()

node = open(sys.argv[3], "w")
node.write("walrus node file\n")
num_nodes = len(v)
node.write("num_nodes " + str(num_nodes) + "\n")
node.write("num_surface_nodes " + str(num_nodes) + "\n")
node.write("num_interior_nodes 0\n")
node.write("surface nodes:\n")
for vert in v:
	node.write(vert)
node.write("interior nodes:\n")
node.close()
