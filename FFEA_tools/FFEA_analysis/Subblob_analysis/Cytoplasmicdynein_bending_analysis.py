import sys, os
import FFEA_node, FFEA_pin, FFEA_trajectory
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

if len(sys.argv) != 8:
	sys.exit("Usage: python Cytoplasmicdynein_bending_analysis.py [Stalkhead subblob (.pin)] [Stalkbase subblob (.pin)] [Axial node 1] [Axial node 2] [INPUT .node] [INPUT trajectory .traj] [FFEA scale]")

# Get args
stalkhead_fname = sys.argv[1]
stalkbase_fname = sys.argv[2]
axial_node = [int(sys.argv[3]), int(sys.argv[4])]
stalknode_fname = sys.argv[5]
traj_fname = sys.argv[6]
scale = float(sys.argv[7])

# Build objects
head_pin = FFEA_pin.FFEA_pin(stalkhead_fname)
base_pin = FFEA_pin.FFEA_pin(stalkbase_fname)
node = FFEA_node.FFEA_node(stalknode_fname)
traj = FFEA_trajectory.FFEA_trajectory(traj_fname)
index = [traj.blob[0][0].define_subblob(head_pin.node_index), traj.blob[0][0].define_subblob(base_pin.node_index)]
subblob = [traj.blob[0][0].subblob[index[0]], traj.blob[0][0].subblob[index[1]]]

head_traj = traj.blob[0][0].get_centroid_trajectory(index[0])
base_traj = traj.blob[0][0].get_centroid_trajectory(index[1])
fluc_traj = [base_traj[i] - head_traj[i] for i in range(len(head_traj))]

# Build co-ordinate system
r0 = base_traj[0] - head_traj[0]

xp = r0
xp *= 1.0 / np.linalg.norm(xp)

zp = traj.blob[0][0].frame[0].pos[axial_node[1]] / scale - traj.blob[0][0].frame[0].pos[axial_node[0]] / scale

yp = np.cross(zp, xp)
yp *= 1.0 / np.linalg.norm(yp)

zp = np.cross(xp, yp)
zp *= 1.0 / np.linalg.norm(zp)

# Build initial projections into new coordinate system
r0xy = r0 - (np.dot(r0, zp)) * zp
r0xy *= 1.0 / np.linalg.norm(r0xy)

r0yz = r0 - (np.dot(r0, xp)) * xp
r0yz *= 1.0 / np.linalg.norm(r0yz)

r0zx = r0 - (np.dot(r0, yp)) * yp
r0zx *= 1.0 / np.linalg.norm(r0zx)

# Get fluctuation distribution
thetaxy = []
thetayz = []
thetazx = []

txymean = 0.0
tyzmean = 0.0
tzxmean = 0.0

txymeanerr = 0.0
tyzmeanerr = 0.0
tzxmeanerr = 0.0

txysd = 0.0
tyzsd = 0.0
tzxsd = 0.0

txysderr = 0.0
tyzsderr = 0.0
tzxsderr = 0.0

for f in fluc_traj:

	# Values
	r = f / scale

	rxy = r - (np.dot(r, zp)) * zp
	rxy *= 1.0 / np.linalg.norm(rxy)

	ryz = r - (np.dot(r, xp)) * xp
	ryz *= 1.0 / np.linalg.norm(ryz)

	rzx = r - (np.dot(r, yp)) * yp
	rzx *= 1.0 / np.linalg.norm(rzx)

	r *= 1.0 / np.linalg.norm(r)

	# Include possibility of negative angles
	if np.dot(np.cross(rxy, r0xy), zp) > 0:
		txy = np.arccos(np.dot(r, rxy)) * 180.0 / np.pi
	else:
		txy = -1 * np.arccos(np.dot(r, rxy)) * 180.0 / np.pi

	if np.dot(np.cross(ryz, r0yz), xp) > 0:
		tyz = np.arccos(np.dot(r, ryz)) * 180.0 / np.pi
	else:
		tyz = -1 * np.arccos(np.dot(r, ryz)) * 180.0 / np.pi

	if np.dot(np.cross(rzx, r0zx), yp) > 0:
		tzx = np.arccos(np.dot(r, rzx)) * 180.0 / np.pi
	else:
		tzx = -1 * np.arccos(np.dot(r, rzx)) * 180.0 / np.pi

	thetaxy.append(txy)
	thetayz.append(tyz)
	thetazx.append(tzx)

	# Moments
	txymean += txy
	tyzmean += tyz
	tzxmean += tzx

	txysd += txy * txy
	tyzsd += tyz * tyz
	tzxsd += tzx * tzx

# Finalise moments
txymean *= 1.0 / len(fluc_traj)
tyzmean *= 1.0 / len(fluc_traj)
tzxmean *= 1.0 / len(fluc_traj)

txysd = np.sqrt((txysd / len(fluc_traj)) - (txymean * txymean))
tyzsd = np.sqrt((tyzsd / len(fluc_traj)) - (tyzmean * tyzmean))
tzxsd = np.sqrt((tzxsd / len(fluc_traj)) - (tzxmean * tzxmean))

txymeanerr = txysd * 1.0 / np.sqrt(len(fluc_traj))
tyzmeanerr = tyzsd * 1.0 / np.sqrt(len(fluc_traj))
tzxmeanerr = tzxsd * 1.0 / np.sqrt(len(fluc_traj))
	
# Move distributions to origin
for i in range(len(thetaxy)):
	thetaxy[i] -= txymean
	thetayz[i] -= tyzmean
	thetazx[i] -= tzxmean

# Make some bins

# Plot the Distributions on a single figure
plt.clf()
plt.cla()
n, bins, patches = plt.hist(thetaxy, 50, normed = 1, facecolor='green', alpha=0.75, label="Stalk Base-Head Angular Distribution")
x = np.linspace(-3 * txysd,3 * txysd,50)
y = mlab.normpdf(x, 0, txysd)

l, = plt.plot(x, y, 'r--', linewidth=1, label=r"Best Fit Normal Distribution" + "\n" + r"$  \mu = %5.2f \pm %5.2f$" % (0.0, txymeanerr) + "\n" + r"$  \sigma = %5.2f$" % (txysd))

plt.xlabel('Angle (degrees)')
plt.ylabel('Probability Density')
plt.title('Dynein Stalk Fluctuations')


handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1], loc=1, fontsize=11, fancybox=True, shadow=True)
#plt.grid(True)

plt.show()
