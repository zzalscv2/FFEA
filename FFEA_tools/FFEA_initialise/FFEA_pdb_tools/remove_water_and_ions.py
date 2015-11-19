import sys, os
import FFEA_pdb

if len(sys.argv) != 3:
	sys.exit("Usage: python convert_dynein.py [INPUT .pdb fname] [OUTPUT .pdb fname]")

# Get args
infname = sys.argv[1]
outfname = sys.argv[2]

# Do stuff
pdb = FFEA_pdb.FFEA_pdb(infname)
pdb.reorder_atoms()
pdb.reorder_residues()
pdb.write_to_file(outfname)
