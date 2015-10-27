import sys, os
from math import floor
import FFEA_trajectory
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 4:
	sys.exit("Usage: python FFEA_test_convergence.py [FFEA .traj fname] [Simulation scale] [num_PCA checks]")

# Get args
script_fname = os.path.abspath(sys.argv[0])
traj_fname = os.path.abspath(sys.argv[1])
scale = float(sys.argv[2])
num_test_points = int(sys.argv[3])
total_num_frames = FFEA_trajectory.get_num_frames(traj_fname)
num_frames = np.array([int(floor((total_num_frames * (i + 1)) / num_test_points)) for i in range(num_test_points)])

# Get current directory and other fnames
scriptdirname = os.path.dirname(script_fname)
trajbasename = os.path.splitext(traj_fname)[0]
convert = scriptdirname + "/FFEA_convert_FFEA_to_atomic.py"
out_fname = trajbasename + ".out"
pdb_top_fname = os.path.splitext(out_fname)[0] + "_frame0.pdb"
pdb_traj_fname = trajbasename + ".pdb"
FFEA_top_fname = os.path.splitext(out_fname)[0] + "_frame0.out"
pcz_fname = trajbasename + ".pcz"
eval_fname = trajbasename + ".evals"

# Run PCA multiple times (this could take ages!)
alleigs = []
for i in range(num_test_points):

	# Get PCA requirements
	os.system("python " + convert + " -traj %s -scale %e -out %s -format pdb -frames %d" % (traj_fname, scale, out_fname, num_frames[i]))
	
	# Run the PCA
	os.system("pyPcazip -i %s -t %s -o %s -e 5 -vvv" % (pdb_traj_fname, pdb_top_fname, pcz_fname))

	# Analyse the PCA (get the eigenvalues)
	os.system("pyPczdump -i %s -l -o %s" % (pcz_fname, eval_fname))

	# Get the eigenvalues into a list
	eigsframe = []
	with open(eval_fname, "r") as fin:
		for line in fin.readlines():
			eigsframe.append(float(line))
	
	alleigs.append(eigsframe)	

# Begin plotting
alleigs = np.array(alleigs)
alleigs = alleigs.T
for i in range(len(alleigs)):
	plt.plot(num_frames, alleigs[i])

plt.xlabel("Num Frames Analysed")
plt.ylabel("Eigenvalue")
plt.title("Eigenvalue Convergence")
plt.show()
figname = trajbasename + "_evalconv.png"
plt.savefig(figname)
print("If pyPca/1.4.0 was not loaded, this program didn't work! Load pyPca and run again.")
os.system("rm " + eval_fname)
