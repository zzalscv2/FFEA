import os, sys

if len(sys.argv) != 2:
	sys.exit("Usage: python make_default_LJmatrix.py [OUTPUT .LJ FILE]")

num_vdw_types = 7

LJmat = open(sys.argv[1], "w")
LJmat.write("ffea vdw forcefield params file\n")
LJmat.write("num_vdw_face_types " + str(num_vdw_types) + "\n")

for i in range(num_vdw_types):
	for j in range(num_vdw_types):
		if i == 0 and j == 0:
			vdw_eps = 1e15
		else:
			vdw_eps = 1e12
		
		vdw_req = 1e-9
		line = "(" + str(vdw_eps) + ", " + str(vdw_req) + ") "
		LJmat.write(line)
	LJmat.write("\n")
LJmat.close()
