import sys, os
import FFEA_meas

if len(sys.argv) != 6:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT .meas file] [num_blobs] [num_nodes] [num_frames_to_read] [frame length]")

# Get args
inmeas = sys.argv[1]
basename = os.path.splitext(os.path.abspath(inmeas))[0]
num_blobs = int(sys.argv[2])
num_nodes = [int(sys.argv[3])]
frames_to_read = int(sys.argv[4])
frame_length = float(sys.argv[5])
kT = 4.11e-21

# Get meas object
meas = FFEA_meas.FFEA_meas(inmeas, num_blobs, frames_to_read, frame_length)

# Plot stuff
meas.plot_energies(basename, num_nodes, kT)
