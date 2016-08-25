import sys, os
import numpy as np

if len(sys.argv) != 3:
	sys.exit("Usage: python FFEA_get_PCA_projections.py [INPUT .pcz file] [num_modes]")

# Get args
infile = sys.argv[1]
base, ext = os.path.splitext(os.path.abspath(infile))
num_modes = int(sys.argv[2])
projfname = base + ".proj"
tempprojfname = "temp.proj"

# Assemble projections into arrays. Must only be the first order nodes! Therefore, must renormalise
print("Calculating Projections: Writing to " + os.path.basename(projfname) + "...")

fout = open(projfname, "w")
for i in range(num_modes):
	proj = []
	print("\tProjection " + str(i) + "...")
	os.system("pyPczdump -i %s -p %d -o %s" % (infile, i, tempprojfname))
	fin = open(tempprojfname, "r")
	lines = fin.readlines()
	for l in lines:
		proj.append(float(l))

	for elem in proj:
		fout.write("%6.3f " % (elem))
	fout.write("\n")
	print("\tdone!")
	fin.close()
fout.close()
print("done!")
os.system("rm " + tempprojfname)
