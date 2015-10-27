import sys, os
import numpy as np
import FFEA_topology

if len(sys.argv) != 4:
	sys.exit("Usage: python FFEA_get_eigensystem.py [INPUT .pcz file] [INPUT Original FFEA .top file] [num_modes]")

# Get args
infile = sys.argv[1]
base, ext = os.path.splitext(os.path.abspath(infile))
top = FFEA_topology.FFEA_topology(sys.argv[2])
lin = list(top.get_linear_nodes())
num_modes = int(sys.argv[3])
eigvalfname = base + ".evals"
tempeigvecfname = "temp.evecs"
eigvecfname = base + ".evecs"

# Send eigenvalues to file (easy)
print("Calculating Eigenvalues: Writing to " + os.path.basename(eigvalfname) + "...")
os.system("pyPczdump -i %s -l -o %s" % (infile, eigvalfname))
print("done!")

# Assemble eigenvectors and eigenvalues into arrays. Must only be the first order nodes! Therefore, must renormalise
fout = open(eigvecfname, "w")
print("Calculating Eigenvectors: Writing to " + os.path.basename(eigvecfname) + "...")

for i in range(num_modes):
	eigvec = []
	print("\tEigenvector " + str(i) + "...")
	os.system("pyPczdump -i %s -e %d -o %s" % (infile, i, tempeigvecfname))
	fin = open(tempeigvecfname, "r")
	lines = fin.readlines()
	num_nodes = len(lines) / 3
	for j in range(num_nodes):
		for k in range(3):
			if j in lin:
				eigvec.append(float(lines[3 * j + k]))			

	eigvec = np.array(eigvec)
	eigvec /= np.linalg.norm(eigvec)

	for elem in eigvec:
		fout.write("%6.3f " % (elem))
	fout.write("\n")
	print("\tdone!")
	fin.close()
fout.close()
print("done!")
os.system("rm " + tempeigvecfname)
