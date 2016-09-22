import sys, os
import numpy as np
import copy

if len(sys.argv) < 4:
	sys.exit("Usage: python " + sys.argv[0] + " [INPUT .evecs file 1] [INPUT .evecs file 2] [OUTPUT fname] OPTIONAL{[EIG SWITCH (x,y)]}")

infile = [sys.argv[1], sys.argv[2]]
out_fname = sys.argv[3]
evec_fname = [[],[]]
evecs = [[],[]]

eig_switch = []
if len(sys.argv) > 4:
	for arg in sys.argv[4:]:
		try:
			eig_switch.append([int(a) for a in arg.split(",")])
		except:
			sys.exit("EIG SWITCH input error. Pairs of integer values separated by columns please")

# Get motion eigenvectors
print("Reading eigenvectors from inputs...")
for i in range(2):
	fin = open(infile[i])
	for line in fin.readlines():
		sline = line.split()
		evecs[i].append(np.array([float(s) for s in sline]))
	fin.close()
	#evecs[i] = np.array(evecs[i])

print("done!")

# In one of these, switch the eigenvectors of necessary
if eig_switch != []:
	for e in eig_switch:
		try:
			e1 = copy.copy(evecs[0][e[0]])
			e2 = copy.copy(evecs[0][e[1]])
		except:
			sys.exit("EIG SWITCH values greater than number of eigenvectors present. Please try again")

		evecs[0][e[0]] = e2
		evecs[0][e[1]] = e1

# Do some dot products!
print("Building the Eigenvectors Dot Product Matrix...")
num_modes = len(evecs[0])
dot_prod = np.array([[0.0 for i in range(num_modes)] for j in range(num_modes)])

for i in range(num_modes):
	for j in range(num_modes):
		dot_prod[i][j] = np.dot(evecs[0][i], evecs[1][j])

print("done!")

# Print to file
print("Writing Matrix to file " + os.path.basename(out_fname) + "...")
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
