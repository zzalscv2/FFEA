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
import MDAnalysis as mda
import numpy as np

if (len(sys.argv) != 6):
	sys.exit("Usage: python FFEA_pdb_rotate.py [INPUT .pdb fname] [OUTPUT .pdb fname] [x angle (degrees)] [y angle (degrees)] [z angle (degrees)]")
	
# Get args
infname = sys.argv[1]
outfname = sys.argv[2]

# These should be degrees
x = np.radians(float(sys.argv[3]))
y = np.radians(float(sys.argv[4]))
z = np.radians(float(sys.argv[5]))

# Build matrix
c = np.cos
s = np.sin

R = [[0.0 for i in range(3)] for j in range(3)]
R[0][0] = c(y) * c(z)
R[0][1] = s(x) * s(y) * c(z) - c(x) * s(z)
R[0][2] = c(x) * s(y) * c(z) + s(x) * s(z)
R[1][0] = c(y) * s(z)
R[1][1] = s(x) * s(y) * s(z) + c(x) * c(z)
R[1][2] = c(x) * s(y) * s(z) - s(x) * c(z)
R[2][0] = -1 * s(y)
R[2][1] = s(x) * c(y)
R[2][2] = c(x) * c(y)

# Get md object
u = mda.Universe(infname)
pos = u.atoms.CA.center_of_geometry()
u.atoms.translate(-1 * pos)
u.atoms.rotate(R)
u.atoms.translate(pos)

# Write to file
u.atoms.write(outfname)
