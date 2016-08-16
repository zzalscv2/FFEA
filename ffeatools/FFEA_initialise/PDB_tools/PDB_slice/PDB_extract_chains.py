import FFEA_pdb
import sys, os

if len(sys.argv) < 4:
	sys.exit("Usage: python PDB_extract_chains.py [INPUT .pdb fname] [OUTPUT .pdb fname] [List of chains to keep]")

# Get args
print("Reading args...")
infname = sys.argv[1]
outfname = sys.argv[2]
safechains = []
for arg in sys.argv[3:]:
	safechains.append(arg.upper())
print("done!")

# Read in PDB
print("Reading pdb...")
pdb = FFEA_pdb.FFEA_pdb(infname)
print("done!")
print("Deleting chains...")
for i in range(pdb.num_atoms - 1, -1, -1):
	if pdb.atom[i].chain not in safechains:
		pdb.atom.pop(i)
		pdb.num_atoms -= 1

print("done!")

# Write to file
pdb.write_to_file(outfname)	
