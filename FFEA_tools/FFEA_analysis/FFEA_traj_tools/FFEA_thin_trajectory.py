import sys, os
from math import ceil
import FFEA_trajectory

if len(sys.argv) != 5:
	sys.exit("Usage python " + os.path.basename(sys.argv[0]) + " [FFEA traj fname] [FFEA output traj fname] [frames to read] [Percentage to keep]")

traj_fname = sys.argv[1]
out_fname = sys.argv[2]
frames_to_read = int(sys.argv[3])
thin_percent = float(sys.argv[4])

if thin_percent < 0 or thin_percent > 100:
	sys.exit("Error. Percentage must be between 0 and 100. You used %f\n" % (thin_percent))

if thin_percent < 1:
	verify = raw_input("Percentage was %f. Did you mean %f (y/n)?" % (thin_percent, thin_percent * 100))
	if verify.lower() == "y":
		thin_percent *= 100


frame_rate = ceil(100 / thin_percent)
traj = FFEA_trajectory.FFEA_trajectory(traj_fname, frame_rate = frame_rate, num_frames_to_read = frames_to_read)
traj.write_to_file(out_fname)
