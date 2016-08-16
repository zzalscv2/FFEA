import sys
import FFEA_pdb
import numpy as np

if len(sys.argv) != 10:
	sys.exit("Usage: python PDB_rotate.py [INPUT .pdb fname] [OUTPUT .pdb fname] [About axis(x,y,z)] [At point (x,y,z)] [Angle (degrees)]")

# Get args
inpdbfname = sys.argv[1]
outpdbfname = sys.argv[2]
axis = np.array([float(sys.argv[i + 3]) for i in range(3)])
point = np.array([float(sys.argv[i + 6]) for i in range(3)])
angle = float(sys.argv[9])

# Build pdb of doom
pdb = FFEA_pdb.FFEA_pdb(inpdbfname)
pdb.rotate(point, axis, angle)
pdb.write_to_file(outpdbfname)


