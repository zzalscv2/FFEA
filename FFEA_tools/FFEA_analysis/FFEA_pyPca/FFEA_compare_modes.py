import sys, os
import numpy as np
if len(sys.argv) != 5:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .pcz file 1] [INPUT .pcz file 2] [OUTPUT fname] [num_modes]")

infile = [sys.argv[1], sys.argv[2]]
out_fname = sys.argv[3]
num_modes = int(sys.argv[4])
evec_fname = [[],[]]
evecs = [[],[]]
dot_prod = np.array([[0.0 for i in range(num_modes)] for j in range(num_modes)])

# Get motion eigenvectors
for i in range(num_modes):
	for j in range(2):
		anfname = "f_" + str(j) + "_evec_" + str(i) + ".evec"
		evec_fname[j].append(anfname)
		os.system("pyPczdump --input " + infile[j] + " --evec " + str(i) + " > " + anfname)
		fin = open(anfname, "r")
		lines = fin.readlines()
		fin.close()
		evec = np.array([0.0 for k in range(len(lines))])
		for k in range(len(lines)):
			evec[k] = float(lines[k].strip())
		evecs[j].append(evec)

# Do some dot products!
for i in range(num_modes):
	for j in range(num_modes):
		dot_prod[i][j] = np.dot(evecs[0][i], evecs[1][j])

# Print to file
outfile = open(out_fname, "w")
outfile.write("Compare pyPca modes\n\nEigen Set 1 x Eigen Set 2\n\n")
for i in range(num_modes):
	for j in range(num_modes):
		outfile.write("%6.3f " % (dot_prod[i][j]))

	outfile.write("\n")

outfile.close()
# Remove crap
for i in range(2):
	for fname in evec_fname[i]:
		os.system("rm " + fname)
