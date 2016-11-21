import sys, os
import MDAnalysis as mda

if (len(sys.argv) != 6):
	sys.exit("Usage: python FFEA_pdb_set_position.py [INPUT .pdb fname] [OUTPUT .pdb fname] [x] [y] [z]")
	
# Get args
infname = sys.argv[1]
outfname = sys.argv[2]
pos = [float(i) for i in sys.argv[3:]]

# Build mda universe
u = mda.Universe(infname)

# Translate it to x, y, z
u.atoms.translate(pos)

# Write it to file
u.atoms.write(outfname)
