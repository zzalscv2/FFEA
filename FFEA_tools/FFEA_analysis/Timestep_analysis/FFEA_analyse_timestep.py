import sys, os, glob
from math import *
import numpy as np
import matplotlib.pyplot as plt
import FFEA_meas

def get_FFEA_meas_fnames(top_dir):

	# Get all directories in top_meas_dir
	if(top_dir[-1] == "/"):
		dir_fnames = glob.glob(top_dir + "*")
	else:
		dir_fnames = glob.glob(top_dir + "/*")

	# Get all fnames and timesteps
	timestep = []
	meas_fnames = []
	for dirfname in dir_fnames:
		timestep.append(float(os.path.basename(dirfname)))
		meas_fnames.append(glob.glob(dirfname + "/*"))

	# Sort timesteps into order
	together = zip(timestep, meas_fnames)
	sorted_together =  sorted(together)

	timestep_sorted = [x[0] for x in sorted_together]
	meas_fnames_sorted = [x[1] for x in sorted_together]
	return meas_fnames_sorted, timestep_sorted

# Main program
for arg in sys.argv:
	if arg == "-h" or arg == "--help":
		print "A python script for finding the maximum allowed timestep for a given system of FFEA blobs. Expects a directory containing multiple directories, each with a single run of a single system at a given timestep. The directory name should be the timestep\n"
		sys.exit("Usage python " + sys.argv[0] + " [FFEA top level measurement directory] [OUTPUT timestep analysis fname]\n")		

if len(sys.argv) != 5:
	sys.exit("Usage python " + sys.argv[0] + " [FFEA top level measurement directory] [OUTPUT timestep analysis fname] [Correct kinetic energy] [Correct potential energy]\n")

# Arguments
meas_dir = sys.argv[1]
out_fname = sys.argv[2]
correct_ke = float(sys.argv[3])
correct_pe = float(sys.argv[4])

if len(out_fname.split(".")) == 2:
	out_basename = out_fname.split(".")[0]
else:
	out_basename = out_fname

# Get all fnames and timesteps in meas_dir
meas_fnames, timesteps = get_FFEA_meas_fnames(meas_dir)
num_timesteps = len(timesteps)
num_blobs = len(meas_fnames[0]) - 1
base_fnames = [""] * num_timesteps
for i in range(num_timesteps):
	for fname in meas_fnames[i]:
		if "_world" in fname:
			base_fnames[i] = fname.split("_world")[0]
			break

# Get average energy for each timestep
ke = [np.array([0.0 for i in range(num_timesteps)]) for j in range(num_blobs)]
pe = [np.array([0.0 for i in range(num_timesteps)]) for j in range(num_blobs)]
cke = np.array([correct_ke for i in range(num_timesteps)])
cpe = np.array([correct_pe for i in range(num_timesteps)])

for i in range(num_timesteps):
	meas = FFEA_meas.FFEA_meas(base_fnames[i], num_blobs, 10000)
	meas.calc_avg_energies()
	for j in range(num_blobs):
		ke[j][i], pe[j][i] = meas.blob[j].get_avg_energies()


# Plot some graphs
figindex = 0
for i in range(num_blobs):
	plt.figure(figindex)
	act, = plt.semilogx(timesteps, ke[i])
	perf, = plt.semilogx(timesteps, cke)
	plt.title("Average Kinetic Energy vs Timestep")
	plt.xlabel("Timestep (s)")
	plt.ylabel(r"$\langle E_k \rangle (J)$")
	plt.legend([act, perf], ["Simualtion Averages", "Ideal Values"], loc = 2, prop={'size':12})
	plt.savefig(out_basename + "_blob" + str(i) + "kevdt.jpg")

	figindex += 1
	
	plt.figure(figindex)
	act, = plt.semilogx(timesteps, pe[i])
	perf, = plt.semilogx(timesteps, cpe)
	plt.title("Average Potential Energy vs Timestep")
	plt.xlabel("Timestep (s)")
	plt.ylabel(r"$\langle E_p \rangle (J)$")
	plt.legend([act, perf], ["Simualtion Averages", "Ideal Values"], loc = 2, prop={'size':12})
	plt.savefig(out_basename + "_blob" + str(i) + "pevdt.jpg")

	figindex += 1
