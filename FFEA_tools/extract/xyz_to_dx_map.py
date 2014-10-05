#!/usr/bin/env python

import numpy as np
import os
import sys
import math

if len(sys.argv) != 7:
	sys.exit("Usage: python make_dx_map.py [INPUT 3 COLUMN X Y Z FILE] [NX] [NY] [NZ] [SCALING] [OUTPUT DX MAP FILENAME]")
infile_fname = sys.argv[1]
xn = int(sys.argv[2])
yn = int(sys.argv[3])
zn = int(sys.argv[4])
scaling = float(sys.argv[5])
outfile_fname = sys.argv[6]

print "infile_fname =", infile_fname
print "xn =", xn
print "yn =", yn
print "zn =", zn
print "scaling =", scaling
print "outfile_fname =", outfile_fname

# load array from input file
pos_list = np.loadtxt(infile_fname, skiprows=0, dtype=float)

# convert trajectory to angstroms
pos_list *= 1e10

# get the minimum x,y,z position (this will be the "origin" of the map)
# and the maximum x,y,z position (this will determine the scaling required for xn, yn, zn bins)
minx = float("Inf")
miny = float("Inf")
minz = float("Inf")
maxx = float("-Inf")
maxy = float("-Inf")
maxz = float("-Inf")


pos_list = np.transpose(pos_list)
minx = np.amin(pos_list[0])
miny = np.amin(pos_list[1])
minz = np.amin(pos_list[2])
maxx = np.amax(pos_list[0])
maxy = np.amax(pos_list[1])
maxz = np.amax(pos_list[2])
pos_list = np.transpose(pos_list)

print "Without margin:"
print "minx =", minx
print "miny =", miny
print "minz =", minz

print "maxx =", maxx
print "maxy =", maxy
print "maxz =", maxz

print "With margin:"
minx -= 5
miny -= 5
minz -= 5
maxx += 5
maxy += 5
maxz += 5
print "minx =", minx
print "miny =", miny
print "minz =", minz

print "maxx =", maxx
print "maxy =", maxy
print "maxz =", maxz

scale_x = (maxx - minx) / (xn)
scale_y = (maxy - miny) / (yn)
scale_z = (maxz - minz) / (zn)
print "scale_x =", scale_x
print "scale_y =", scale_y
print "scale_z =", scale_z

# place points in bins of 3d histogram
histo = np.zeros((xn, yn, zn), dtype=float)
num_atoms = 0
for i in range(len(pos_list)):
	atomx = pos_list[i][0]
	atomy = pos_list[i][1]
	atomz = pos_list[i][2]
	binx = math.floor(((atomx - minx) / scale_x))
	biny = math.floor(((atomy - miny) / scale_y))
	binz = math.floor(((atomz - minz) / scale_z))
	histo[binx][biny][binz] += 1
	num_atoms += 1

# normalise map and scale it
histo = (histo / num_atoms) * scaling

dx = open(outfile_fname, "w")
dx.write("# Wassup bruv?\n")
dx.write("object 1 class gridpositions counts " + str(xn) + " " + str(yn) + " " + str(zn) + "\n")
dx.write("origin " + str(minx) + " " + str(miny) + " " + str(minz) + "\n")
dx.write("delta " + str(scale_x) + " 0.0 0.0\n")
dx.write("delta 0.0 " + str(scale_y) + " 0.0\n")
dx.write("delta 0.0 0.0 " + str(scale_z) + "\n")
dx.write("object 2 class gridconnections counts " + str(xn) + " " + str(yn) + " " + str(zn) + "\n")
dx.write("object 3 class array type double rank 0 times " + str(xn * yn * zn) + " data follows\n")
i = 0
for x in range(xn):
	for y in range(yn):
		for z in range(zn):
			dx.write(str(histo[x][y][z]))
			i += 1
			if i == 3:
				i = 0
				dx.write("\n")
			else:
				dx.write(" ")
if i != 3:
	dx.write("\n")
dx.write("attribute \"dep\" string \"positions\"\n")
dx.write("object \"regular positions regular connections\" class field\n")
dx.write("component \"positions\" value 1\n")
dx.write("component \"connections\" value 2\n")
dx.write("component \"data\" value 3\n")
dx.close()
