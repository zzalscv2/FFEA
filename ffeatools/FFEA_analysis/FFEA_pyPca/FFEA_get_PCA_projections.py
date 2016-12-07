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
