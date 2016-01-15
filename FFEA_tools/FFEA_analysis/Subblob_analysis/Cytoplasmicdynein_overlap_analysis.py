import sys, os
import numpy as np
import FFEA_node, FFEA_pin, FFEA_trajectory
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

if len(sys.argv) != :
	sys.exit("Usage: python Cytoplasmicdynein_overlap_analysis.py [INPUT .traj] [Subblob 1 (.pin)] [Subblob 2 (.pin)]")

# Get args
trajfname = sys.argv[1]
pinfname = [sys.argv[2], sys.argv[3]]

# Get objects
traj = FFEA_trajectory.FFEA_trajectory(trajfname)
pin = [FFEA_pin.FFEA_pin(p) for p in pinfname]

# Build subblobs
traj.blob[0][0].define_subblob(pin[0].node_index)
traj.blob[1][0].define_subblob(pin[1].node_index)
motor = [traj.blob[0][0].subblob[0], traj.blob[1][0].subblob[0]]

# Get trajectories of motors only!
motortraj = np.array([traj.blob[0][0].get_centroid_trajectory(0), traj.blob[0][0].get_centroid_trajectory(0)])
num_frames = len(motortraj[0])

# Now, get the distribution of separations in the axial (z) direction
r = motortraj[1] - motortraj[0]
print motortraj[1]
print motortraj[0]
print r[0]
sys.exit()
zsep = np.array([np.dot((r[i] - r[0]), np.array([0.0,0.0,1.0])) * 1e9 for i in range(num_frames)])
zsepmean = 0.0
zsepsd = 0.0
for sep in zsep:
	zsepmean += sep
	zsepsd += sep * sep

zsepmean *= 1.0 / num_frames
zsepsd = np.sqrt(zsepsd * 1.0 / num_frames) - zsepmean)
zsepmeanerr = zsepsd * 1.0 / sqrt(num_frames)

# Plot the distribution and fit a histogram
n, bins, patches = plt.hist(zsep, 10, normed=1, facecolor='g', alpha=0.75, label="Seperation")
y = mlab.normpdf(zsep, zsepmean, zsepsd)
plt.plot(zsep, y, 'r--', linewidth=1, label=r"Best Fit Normal Distribution" + "\n" + r"$  \mu = %5.2f \pm %5.2f$" % (0.0, zsepmeanerr) + "\n" + r"$  \sigma = %5.2f$" % (zsepsdsd))

# Fix graph properties
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], loc=1, fontsize=11, fancybox=True, shadow=True)
plt.xlabel("Distance (nm)")
plt.ylabel("Probability")
plt.title("Cytoplasmic Dynein Motor Domain Axial Separation")
