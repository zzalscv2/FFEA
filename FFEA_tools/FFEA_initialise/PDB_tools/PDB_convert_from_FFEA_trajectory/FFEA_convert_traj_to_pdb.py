import sys, os
import FFEA_trajectory, FFEA_pdb

if len(sys.argv) != 5:
	sys.exit("Usage: python FFEA_convert_traj_to_pdb.py [INPUT FFEA traj fname (.out)] [OUTPUT .pdb fname] [num_frames_to_read] [FFEA scale]")

# Get args
infname = sys.argv[1]
outfname = sys.argv[2]
num_frames_to_read = int(sys.argv[3])
scale = 1.0 / float(sys.argv[4])	# Invert as pdb is in angstroms

# Build objects
traj = FFEA_trajectory.FFEA_trajectory(infname)
pdb = FFEA_pdb.FFEA_pdb("")
pdb.build_from_traj(traj, scale = scale)
pdb.write_to_file(outfname)
