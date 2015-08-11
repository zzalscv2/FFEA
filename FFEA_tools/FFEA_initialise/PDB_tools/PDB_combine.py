import sys
import FFEA_pdb
import numpy as np

if len(sys.argv) != 4:
	sys.exit("Usage: python PDB_rotate.py [INPUT .pdb fname 1] [INPUT .pdb fname 2] [OUTPUT .pdb fname] ")

# Get args
inpdbfname = [sys.argv[1], sys.argv[2]]
outpdbfname = sys.argv[3]

# Build pdb of doom
pdb = [FFEA_pdb.FFEA_pdb(fname) for fname in inpdbfname]
pdb[0].add_pdb(pdb[1])
pdb[0].write_to_file(outpdbfname)

