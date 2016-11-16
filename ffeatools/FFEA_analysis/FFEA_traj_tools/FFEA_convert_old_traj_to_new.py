import sys, os
import FFEA_trajectory

if len(sys.argv) != 3:
	sys.exit("python convert_old_traj_to_new.py [FFEA old traj fname] [FFEA new traj fname]")

# Make a traj
traj = FFEA_trajectory.FFEA_trajectory(sys.argv[1])
traj.write_to_file(sys.argv[2])
