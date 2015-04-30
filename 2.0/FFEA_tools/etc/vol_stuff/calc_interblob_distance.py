import sys, math;
from math import *
if len(sys.argv) != 3:
	sys.exit("Usage: calc_interblob_distance.py [INPUT .VOL FILE 1] [INPUT .VOL FILE 2]")

inputvol1 = open(sys.argv[1], "r")
inputvol2 = open(sys.argv[2], "r")

lines1 = inputvol1.readlines()
lines2 = inputvol2.readlines()

lines1index = lines1.index("points\n")
lines2index = lines2.index("points\n")

num_nodes1 = int(lines1[lines1index + 1])
num_nodes2 = int(lines2[lines2index + 1])

for i in range(len(lines1) - 1, lines1index + num_nodes1, -1):
	lines1.pop(i)

for i in range(lines1index + 1, -1, -1):
	lines1.pop(i)

for i in range(len(lines2) - 1, lines2index + num_nodes2, -1):
	lines2.pop(i)

for i in range(lines2index + 1, -1, -1):
	lines2.pop(i)

#Calculating Centroids
centroid1_x = 0.0
centroid1_y = 0.0
centroid1_z = 0.0
centroid2_x = 0.0
centroid2_y = 0.0
centroid2_z = 0.0
for line in lines1:
	sline = line.split()
	centroid1_x = centroid1_x + float(sline[0])
	centroid1_y = centroid1_y + float(sline[1])
	centroid1_z = centroid1_z + float(sline[2])

centroid1_x = centroid1_x / len(lines1)
centroid1_y = centroid1_y / len(lines1)
centroid1_z = centroid1_z / len(lines1)

for line in lines2:
	sline = line.split()
	centroid2_x = centroid2_x + float(sline[0])
	centroid2_y = centroid2_y + float(sline[1])
	centroid2_z = centroid2_z + float(sline[2])

centroid2_x = centroid2_x / len(lines1)
centroid2_y = centroid2_y / len(lines1)
centroid2_z = centroid2_z / len(lines1)

print "Centroid Separation = (" + str(centroid2_x - centroid1_x) + " , " + str(centroid2_y - centroid1_y) + " , " + str(centroid2_z - centroid1_z) + ")\n"
print "\t\t= " + str(sqrt(pow(centroid2_x - centroid1_x, 2) + pow(centroid2_y - centroid1_y, 2) + pow(centroid2_z - centroid1_z, 2)))

#Calculating Min Separation
distance_mag_old = float("inf")
for line1 in lines1:
	sline1 = line1.split()
	for line2 in lines2:
		sline2 = line2.split()
		distance_mag_new = sqrt(pow(float(sline2[0]) - float(sline1[0]), 2) + pow(float(sline2[1]) - float(sline1[1]), 2) + pow(float(sline2[2]) - float(sline1[2]), 2))
		if distance_mag_new < distance_mag_old:
			distance_mag_old = distance_mag_new

print "Minimum Separation = " + str(distance_mag_old) + "\n"








