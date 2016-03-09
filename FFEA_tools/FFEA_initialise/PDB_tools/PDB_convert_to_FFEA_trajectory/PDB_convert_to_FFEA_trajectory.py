import sys, os
import FFEA_pdb, FFEA_trajectory

if len(sys.argv) != 4:
	print("Usage: python PDB_convert_to_FFEA_trajectory.py [INPUT .pdb fname] [OUTPUT FFEA traj fname] [Simulation Scale]")

# Get args
pdbfname = sys.argv[1]
outfname = sys.argv[2]
scale = float(sys.argv[3])

# Build objects
pdb = FFEA_pdb.FFEA_pdb(pdbfname)
traj = FFEA_trajectory.FFEA_trajectory("", load_all = 0)

# Write to FFEA_traj
traj.build_from_pdb(pdb, scale = scale)
traj.write_to_file(outfname)
