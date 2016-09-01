import sys, os
from FFEA_trajectory import get_num_frames as gnf

if len(sys.argv) != 2:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [FFEA traj .ftj]")

print gnf(sys.argv[1])

