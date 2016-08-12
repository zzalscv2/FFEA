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
