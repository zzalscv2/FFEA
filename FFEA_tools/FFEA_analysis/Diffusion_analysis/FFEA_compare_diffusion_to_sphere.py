import sys, os
from math import *
import numpy as np
import matplotlib.pyplot as plt

# Main Program
for arg in sys.argv:
	if arg == "-h" or arg == "--help":
		print "A python script for comparing the diffusion of FFEA blobs to equivalent spheres. Expects a FFEA diffusion analysis file (from 'FFEA_analyse_diffusion') and some sphere properties.\n"
		sys.exit("Usage: python " + sys.argv[0]  + " [FFEA diffusion analysis file (.diffsion)] [sphere radius (m)] [external viscosity (m^2s^-1)] [kbT (J)]")

if len(sys.argv) < 6:
	sys.exit("Usage: python " + sys.argv[0] + " [FFEA diffusion analysis file (.diffsion)] [sphere radius (m)] [sphere_density (kg/m^3)] [external viscosity (m^2s^-1)] [kbT (J)] {[num_nodes]}")
		
infname = sys.argv[1]
out_basename = sys.argv[1].split(".")[0]
radius = float(sys.argv[2])
density = float(sys.argv[3])
mass = (4.0 / 3.0) * np.pi * pow(radius, 3) * density
viscosity = float(sys.argv[4])
kbT = float(sys.argv[5])

# Read diffusion file
num_nodes = []
fin = open(infname, "r")
fin.readline()

num_blobs = int(fin.readline().split(" ")[2])
sline = fin.readline().split(" ")
for i in range(num_blobs):
	num_nodes.append(int(sys.argv[6 + i]))

num_frames = int(fin.readline().split(" ")[2])

time = np.array([0.0 for i in range(num_frames - 1)])
x2 = np.array([[[0.0 for i in range(num_frames - 1)] for j in range(4)] for k in range(num_blobs)])
lnx2 = np.array([[[0.0 for i in range(num_frames - 1)] for j in range(4)] for k in range(num_blobs)])
x2err = np.array([[[0.0 for i in range(num_frames - 1)] for j in range(4)] for k in range(num_blobs)])
lnx2err = np.array([[[0.0 for i in range(num_frames - 1)] for j in range(4)] for k in range(num_blobs)])

fin.readline()
fin.readline()
fin.readline()
for i in range(1, num_frames):

	sline = fin.readline().split()
	time[i - 1] = float(sline[0])

	for j in range(num_blobs):
		for k in range(4):
			x2[j][k][i - 1] = sline[8 * j + 2 * k + 1]
			lnx2[j][k][i - 1] = np.log(x2[j][k][i - 1])
			x2err[j][k][i - 1] = sline[8 * j + 2 * k + 1 + 1]
			lnx2err[j][k][i - 1] = x2err[j][k][i - 1] / x2[j][k][i - 1]


# Plot some graphs
#for i in range(num_blobs):
#	diff_const = 2 * kbT / (6 * np.pi * viscosity * radius)
#	perfect = np.array([t * diff_const for t in time])
#	perfectr = np.array([3 * t * diff_const for t in time])
#	for j in range(4):
#		plt.figure(4 * i + j)
#		plt.errorbar(time, x2[i][j])
#		if j == 3:
#			plt.plot(time, perfectr)
#		else:
#			plt.plot(time, perfect)
#
#		plt.title("Blob %d Diffusion in $x_%d$" % (i,j))
#		plt.xlabel("Time (s)")
#		plt.ylabel(r"$\langle x_%d ^2\rangle (m^2)$" % (j))
#		plt.savefig(out_basename + "_blob" + str(i) + "_x" + str(j) + "^2.jpg")

# Plot some graphs (new)
for i in range(num_blobs):
	diff_const = 2 * kbT / (6 * np.pi * viscosity * radius)
	perfectmass = np.array([2 * np.log(t) + np.log(kbT / mass) for t in time])
	perfectmassr = np.array([2 * np.log(t) + np.log(3 * kbT / mass) for t in time])
	perfectvisc = np.array([np.log(t) + np.log(diff_const) for t in time])
	perfectviscr = np.array([np.log(t) + np.log(3 * diff_const) for t in time])
	for j in range(4):
		plt.figure(4 * i + j)
		dh = plt.errorbar(np.log(time), lnx2[i][j])
		if j == 3:
			pmh, = plt.plot(np.log(time), perfectmassr)
			pvh, = plt.plot(np.log(time), perfectviscr)
		else:
			pmh, = plt.plot(np.log(time), perfectmass)
			pvh, = plt.plot(np.log(time), perfectvisc)

		plt.legend([dh, pmh, pvh], ["Simulation Diffusion", "Ideal Sphere Ballistic Diffusion", "Ideal Sphere Viscous Diffusion"], loc = 2, prop={'size':10})
		plt.xlabel("ln(time)")
		if j == 0:
			plt.title("Blob %d Diffusion in x" % (i))
			plt.ylabel(r"ln($\langle x ^2\rangle$)")
			plt.savefig(out_basename + "_blob" + str(i) + "_x^2.jpg")
		elif j == 1:
			plt.title("Blob %d Diffusion in y" % (i))
			plt.ylabel(r"ln($\langle y ^2\rangle$)")
			plt.savefig(out_basename + "_blob" + str(i) + "_y^2.jpg")
		elif j == 2:
			plt.title("Blob %d Diffusion in z" % (i))
			plt.ylabel(r"ln($\langle z ^2\rangle$)")
			plt.savefig(out_basename + "_blob" + str(i) + "_z^2.jpg")
		elif j == 3:
			plt.title("Blob %d Diffusion in r" % (i))
			plt.ylabel(r"ln($\langle r ^2\rangle$)")
			plt.savefig(out_basename + "_blob" + str(i) + "_r^2.jpg")
