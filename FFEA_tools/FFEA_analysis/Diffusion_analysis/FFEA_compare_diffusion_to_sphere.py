import sys, os
from math import *
import numpy as np
import matplotlib.pyplot as plt

# Main Program
for arg in sys.argv:
	if arg == "-h" or arg == "--help":
		print "A python script for comparing the diffusion of FFEA blobs to equivalent spheres. Expects a FFEA diffusion analysis file (from 'FFEA_analyse_diffusion') and some sphere properties.\n"
		sys.exit("Usage: python " + sys.argv[0]  + " [FFEA diffusion analysis file (.diffsion)] [sphere radius (m)] [external viscosity (m^2s^-1)] [kbT (J)]")

if len(sys.argv) != 5:
	sys.exit("Usage: python " + sys.argv[0] + " [FFEA diffusion analysis file (.diffsion)] [sphere radius (m)] [external viscosity (m^2s^-1)] [kbT (J)]")
		
infname = sys.argv[1]
out_basename = sys.argv[1].split(".")[0]
radius = float(sys.argv[2])
viscosity = float(sys.argv[3])
kbT = float(sys.argv[4])

# Read diffusion file
fin = open(infname, "r")
fin.readline()
sline = fin.readline().split(" ")
num_blobs = int(sline[2])
num_frames = int(sline[5])
time = np.array([0.0 for i in range(num_frames)])
x2 = np.array([[[0.0 for i in range(num_frames)] for j in range(4)] for k in range(num_blobs)])
x2err = np.array([[[0.0 for i in range(num_frames)] for j in range(4)] for k in range(num_blobs)])
fin.readline()
fin.readline()
for i in range(num_frames):
	sline = fin.readline().split()
	time[i] = float(sline[0])
	for j in range(num_blobs):
		for k in range(4):
			x2[j][k][i] = sline[8 * j + 2 * k + 1]
			x2err[j][k][i] = sline[8 * j + 2 * k + 1 + 1]


# Plot some graphs
diff_const = kbT / (6 * np.pi * viscosity * radius)
perfect = np.array([t * diff_const for t in time])
perfectr = np.array([3 * t * diff_const for t in time])
for i in range(num_blobs):
	for j in range(4):
		plt.figure(4 * i + j)
		plt.errorbar(time, x2[i][j])
		if j == 3:
			plt.plot(time, perfectr)
		else:
			plt.plot(time, perfect)

		plt.title("Blob %d Diffusion in $x_%d$" % (i,j))
		plt.xlabel("Time (s)")
		plt.ylabel(r"$\langle x_%d ^2\rangle (m^2)$" % (j))
		plt.savefig(out_basename + "_blob" + str(i) + "_x" + str(j) + "^2.jpg")
