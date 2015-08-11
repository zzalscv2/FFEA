import sys
import FFEA_pdb

if len(sys.argv) != 3:
	sys.exit("Usage: python PDB_renumber.py [INPUT .pdb] [OUTPUT .pdb]")

inpdb = FFEA_pdb.FFEA_pdb(sys.argv[1])
inpdb.write_to_file(sys.argv[2])
