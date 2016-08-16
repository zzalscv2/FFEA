import sys, os
from math import ceil
import FFEA_traj

if len(sys.argv) != 6:
	sys.exit("Usage python " + os.path.basename(sys.argv[0]) + " [FFEA traj fname] [FFEA output traj fname] [frames to read] [First frame] [Last frame]")

traj_fname = sys.argv[1]
out_fname = sys.argv[2]
frames_to_read = int(sys.argv[3])
first_frame = int(sys.argv[4])
last_frame = int(sys.argv[5])

traj = FFEA_traj.FFEA_traj(traj_fname, frames_to_read, first_frame, last_frame, 1)
traj.write_traj_to_file(out_fname)
