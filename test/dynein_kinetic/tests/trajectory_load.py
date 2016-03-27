import sys, os, time
import FFEA_trajectory

if len(sys.argv) != 2 and len(sys.argv) != 4:
	sys.exit("Usage: python traj_load_speed.py [INPUT FFEA traj (.out)] [INPUT FFEA surf 1] [INPUT FFEA surf 2]")

# Get args
fname = sys.argv[1]

# Get beginnings of traj
traj = FFEA_trajectory.FFEA_trajectory(fname, load_all = 0)

# Now time the loading of all frames
start = time.clock()
i = 0
while traj.load_frame() != 0:
	i += 1
	if i % 10 == 0:
		print i, " frames read"	
	
	if i == 1000:
		break

print "Total time to read = ", time.clock() - start, "s"


