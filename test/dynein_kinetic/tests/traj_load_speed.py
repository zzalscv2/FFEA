import sys, os, time
import FFEA_trajectory

if len(sys.argv) != 2:
	sys.exit("Usage: python traj_load_speed.py [INPUT FFEA traj (.out)]")

# Get args
fname = sys.argv[1]

# Get beginnings of traj
traj = FFEA_trajectory.FFEA_trajectory(fname, load_all = 0)

# Now time the loading of all frames
start = time.clock()

