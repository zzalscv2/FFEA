import sys, os, glob
from math import *
import numpy as np
import matplotlib.pyplot as plt
import FFEA_traj

# Extract FFEA traj fnames using glob
def get_FFEA_traj_fnames(dir_fname):

	# Get all files in traj_dir
	if(dir_fname[-1] == "/"):
		traj_fnames = glob.glob(traj_dir + "*")
	else:
		traj_fnames = glob.glob(traj_dir + "/*")

	# If empty, no files
	if len(traj_fnames) == 0:
		sys.exit("Error. No files contained within the directory '" + traj_dir + "'\n")

	# Check they are all FFEA_traj_files. Delete those that aren't
	print "Removing files that aren't FFEA trajectories..."
	remove_list = []
	for fname in traj_fnames:
		try:
			fin = open(fname, "r")

		except IOError:
			remove_list.append(fname)
			continue

		if fin.readline().strip() != "FFEA_trajectory_file":
			remove_list.append(fname)

		fin.close()

	traj_fnames = [fname for fname in traj_fnames if fname not in remove_list]
	remove_list = None

	print "done. " + str(len(traj_fnames)) + " files remaining."
	return traj_fnames

# Main program
for arg in sys.argv:
	if arg == "-h" or arg == "--help":
		print "A python script for analysing diffusion of FFEA blobs. Expects a directory containing multiple runs of a single system.\n"
		sys.exit("Usage python " + sys.argv[0] + " [FFEA trajectory directory]\n")		

if len(sys.argv) != 5:
	sys.exit("Usage python " + sys.argv[0] + " [FFEA trajectory directory] [OUTPUT diffusion analysis fname] [Time per frame (s)] [Print graphs? (y/n)]\n")

# Arguments
traj_dir = sys.argv[1]
out_fname = sys.argv[2]
if len(out_fname.split(".")) == 2:
	out_basename = out_fname.split(".")[0]
else:
	out_basename = out_fname
	out_fname = out_basename + ".diffusion"

time_per_frame = float(sys.argv[3])
print_graphs = sys.argv[4]

# Get all files in traj_dir
traj_fnames = get_FFEA_traj_fnames(traj_dir)
num_trajs = len(traj_fnames)

# Get a centroid trajectory for all trajectories
centroid_traj = []
for fname in traj_fnames:
	centroid_traj.append(FFEA_traj.FFEA_centroid_traj(FFEA_traj.FFEA_traj(fname, 1000000)))
num_blobs = centroid_traj[0].num_blobs
num_frames = centroid_traj[0].num_frames

# Use centroid trajs to calculate <x^2> vs frame number
x2 = np.zeros_like(centroid_traj[0].pos)
x2err = np.zeros_like(centroid_traj[0].pos)
for i in range(num_trajs):
	for j in range(num_blobs):
		x2[j] += np.power(centroid_traj[i].pos[j], 2)
		x2err[j] += np.power(centroid_traj[i].pos[j], 4)

for i in range(num_blobs):
	x2[i] /= num_trajs
	x2err[i] /= num_trajs

	# Set initial x2err to x2 to avoid rounding errors in sqrt
	for j in range(4):
		x2err[i][j][0] = np.power(x2[i][j][0], 2)

	x2err[i] = np.sqrt((x2err[i] - np.power(x2[i], 2)) / num_trajs) 

# Write to file
fout = open(out_fname, "w")
fout.write("FFEA Diffusion analysis\nnum_blobs = " + str(num_blobs) + " num_frames = " + str(num_frames) + "\n\nTime (s)\t")
for i in range(num_blobs):
	fout.write("Blob " + str(i) + "(x2+-dx2  y2+-dy2  z2+-dz2  r2+-dr2)\t")
fout.write("\n")
for i in range(num_frames):
	fout.write("%5.6e\t" % (i * time_per_frame))
	for j in range(num_blobs):
		for k in range(4):
			fout.write("%5.2e %5.2e " % (x2[j][k][i], x2err[j][k][i]))
		fout.write("\t")
	fout.write("\n")

fout.close()

# Plot stuff too maybe
if print_graphs == "y" or print_graphs == "Y":
	frame = np.array([i * time_per_frame for i in range(num_frames)])
	frameerr = np.zeros_like(frame)
	for i in range(num_blobs):
		for j in range(4):
			plt.figure(4 * i + j)
			plt.errorbar(frame, x2[i][j])
			plt.title("Blob %d Diffusion in $x_%d$" % (i,j))
			plt.xlabel("Time (s)")
			plt.ylabel(r"$\langle x_%d ^2\rangle (m^2)$" % (j))
			plt.savefig(out_basename + "_blob" + str(i) + "_x" + str(j) + "^2.jpg")
