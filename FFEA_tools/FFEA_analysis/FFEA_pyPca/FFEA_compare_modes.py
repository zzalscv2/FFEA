import sys, os
import numpy as np
if len(sys.argv) != 4:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .evecs file 1] [INPUT .evecs file 2] [OUTPUT fname]")

infile = [sys.argv[1], sys.argv[2]]
out_fname = sys.argv[3]
evec_fname = [[],[]]
evecs = [[],[]]

# Get motion eigenvectors
print("Reading eigenvectors from inputs...")
for i in range(2):
	fin = open(infile[i])
	for line in fin.readlines():
		sline = line.split()
		evecs[i].append([float(s) for s in sline])
	fin.close()
	evecs[i] = np.array(evecs[i])

print("done!")

# Do some dot products!
print("Building the Eigenvectors Dot Product Matrix...")
num_modes = len(evecs[0])
dot_prod = np.array([[0.0 for i in range(num_modes)] for j in range(num_modes)])

for i in range(num_modes):
	for j in range(num_modes):
		dot_prod[i][j] = np.dot(evecs[0][i], evecs[1][j])

print("done!")

# Print to file
print("Writing Matrix to file " + os.path.basename(outfile) + "...")
outfile = open(out_fname, "w")
outfile.write("Compare pyPca modes\n\nEigen Set 1 x Eigen Set 2\n\n")
for i in range(num_modes):
	for j in range(num_modes):
		outfile.write("%6.3f " % (dot_prod[i][j]))

	outfile.write("\n")

outfile.close()
print("done!")

# Remove crap
for i in range(2):
	for fname in evec_fname[i]:
		os.system("rm " + fname)

# Plot the heatmap thing
print("Visualising as a heatmap...")
os.system(os.path.dirname("python " + os.path.abspath(sys.argv[0])) + "/plot_eigensystem_comparison.py " + out_fname)
print("done!")
