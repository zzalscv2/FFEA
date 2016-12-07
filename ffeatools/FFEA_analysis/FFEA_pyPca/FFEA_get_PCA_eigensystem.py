# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

import sys, os
import numpy as np
import FFEA_topology

if len(sys.argv) < 4:
	sys.exit("Usage: python FFEA_get_eigensystem.py [INPUT .pcz file] [INPUT .top file] [num_modes] OPTIONAL {[OUTPUT fname]}")

# Get args
infile = sys.argv[1]
if len(sys.argv) == 5:
	base, ext = os.path.splitext(os.path.abspath(sys.argv[4]))
else:
	base, ext = os.path.splitext(os.path.abspath(infile))

top = FFEA_topology.FFEA_topology(sys.argv[2])
num_modes = int(sys.argv[3])
eigvalfname = base + ".evals"
tempeigvecfname = base + "_temp.evecs"
eigvecfname = base + ".evecs"

lin = top.get_linear_nodes()

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
