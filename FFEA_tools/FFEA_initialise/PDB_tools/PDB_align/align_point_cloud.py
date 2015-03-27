import sys, os
import math
import random

def read_pdb_frame(file):
        while True:
                line = file.readline()
                if line == "":
                        break
                line = line.replace("-", " -")
                line = line.split()
                if "REMARK" in line[0]:
                        continue
                if "MODEL" in line[0]:
                        continue
                if "ATOM" in line[0]:
                        frame = []
                        while True:
                                if "ENDMDL" in line[0]:
                                        break
                                if "TER" in line[0]:
                                        break
                                if "END" in line[0]:
                                        break
                                frame.append([float(line[5]), float(line[6]), float(line[7])])
                                line = (file.readline())
                                line = line.replace("-", " -")
                                line = line.split()
                        return len(frame), frame
        return -1, None

def get_centroid(num_atoms, frame):
	centroid = [0, 0, 0]
	for i in xrange(num_atoms):
		centroid[0] += frame[i][0]
		centroid[1] += frame[i][1]
		centroid[2] += frame[i][2]
	centroid[0] /= num_atoms
	centroid[1] /= num_atoms
	centroid[2] /= num_atoms
	return centroid

def construct_rotation_matrix(alpha, beta, gamma):
	rotation =	[\
				[math.cos(alpha) * math.cos(beta), math.cos(alpha) * math.sin(beta) * math.sin(gamma) - math.sin(alpha) * math.cos(gamma), math.cos(alpha) * math.sin(beta) * math.cos(gamma) + math.sin(alpha) * math.sin(gamma)],\
				[math.sin(alpha) * math.cos(beta), math.sin(alpha) * math.sin(beta) * math.sin(gamma) + math.cos(alpha) * math.cos(gamma), math.sin(alpha) * math.sin(beta) * math.cos(gamma) - math.cos(alpha) * math.sin(gamma)],\
				[-math.sin(beta), math.cos(beta) * math.sin(gamma), math.cos(beta) * math.cos(gamma)]\
			]
	return rotation

def mat_mult(M, v):
        result = [0, 0, 0]
        for i in xrange(3):
                for j in xrange(3):
                        result[i] += M[i][j]*v[j]
        return result

def get_diff_penalty(frame_large, frame_small, dx, dy, dz, centroid_large, centroid_small, yaw, pitch, roll):
	rotation_matrix = construct_rotation_matrix(yaw, pitch, roll)
	E = 0
	for p in frame_small:
		new_p = mat_mult(rotation_matrix, p)
		new_p[0] += centroid_small[0] + dx
		new_p[1] += centroid_small[1] + dy
		new_p[2] += centroid_small[2] + dz
		for q in frame_large:
			E += math.sqrt((new_p[0] - (q[0] + centroid_large[0]))**2 + (new_p[1] - (q[1] + centroid_large[1]))**2 + (new_p[2] - (q[2] + centroid_large[2]))**2)
	return E

if len(sys.argv) != 4:
	sys.exit("Usage: python rotate_pdb.py [INPUT LARGE POINT CLOUD PDB] [INPUT SMALL POINT CLOUD PDB] [OUTPUT SMALL POINT CLOUD PDB]")

MAX_TRANSLATION_STEP = 1
MAX_ROTATION_STEP = math.pi/90.0

pdb_large = open(sys.argv[1], "r")
num_atoms_large, frame_large = read_pdb_frame(pdb_large)
pdb_large.close()

pdb_small = open(sys.argv[2], "r")
num_atoms_small, frame_small = read_pdb_frame(pdb_small)
pdb_small.close()

centroid_large = get_centroid(num_atoms_large, frame_large)
print "centroid_large =", centroid_large

centroid_small = get_centroid(num_atoms_small, frame_small)
print "centroid_small =", centroid_small

# start off by translating the small cloud into the large cloud such that they have the same centroid
dx = centroid_large[0] - centroid_small[0]
dy = centroid_large[1] - centroid_small[1]
dz = centroid_large[2] - centroid_small[2]
print "Initial guess at translation:"
print "dx =", dx
print "dy =", dy
print "dz =", dz

# initial guess at rotation is no rotation
yaw = 0
pitch = 0
roll = 0

# centre the clouds on the origin
frame_large = [[x[0] - centroid_large[0], x[1] - centroid_large[1], x[2] - centroid_large[2]] for x in frame_large]
frame_small = [[x[0] - centroid_small[0], x[1] - centroid_small[1], x[2] - centroid_small[2]] for x in frame_small]

# get starting "energy" of the configuration
E = get_diff_penalty(frame_large, frame_small, dx, dy, dz, centroid_large, centroid_small, yaw, pitch, roll)
print "Difference penalty E =", E

T = E
max_T = T
for i in xrange(3000):
	print i, "/ 3000"
	# Make a random translation and a random rotation
	ddx = (random.random() - .5) * MAX_TRANSLATION_STEP
	ddy = (random.random() - .5) * MAX_TRANSLATION_STEP
	ddz = (random.random() - .5) * MAX_TRANSLATION_STEP
	dyaw = (random.random() - .5) * MAX_ROTATION_STEP
	dpitch = (random.random() - .5) * MAX_ROTATION_STEP
	droll = (random.random() - .5) * MAX_ROTATION_STEP

	new_E = get_diff_penalty(frame_large, frame_small, dx + ddx, dy + ddy, dz + ddz, centroid_large, centroid_small, yaw + dyaw, pitch + dpitch, roll + droll)
	dE = new_E - E
#	print new_E, dE
#	p = math.exp(-dE/T)
#	if random.random() < p:
	if dE < 0:
		E = new_E
		dx += ddx
		dy += ddy
		dz += ddz
		yaw += dyaw
		pitch += dpitch
		roll += droll
		T -= max_T/20000
		print "T =", T, "E =", E
		print dx, dy, dz, yaw, pitch, roll
#		if T < .1 * max_T:
#			print "T has cooled down to 10%. Ending early."
#			break

pdb_out = open(sys.argv[3], "w")
i = 0
rotation_matrix = construct_rotation_matrix(yaw, pitch, roll)
for pos in frame_small:
	new_pos = mat_mult(rotation_matrix, pos)
	new_pos[0] += centroid_small[0] + dx
	new_pos[1] += centroid_small[1] + dy
	new_pos[2] += centroid_small[2] + dz

	ridiculous_pos_x = ("%.3f" % (new_pos[0])).rjust(12, " ")
	ridiculous_pos_y = ("%.3f" % (new_pos[1])).rjust(8, " ")
	ridiculous_pos_z = ("%.3f" % (new_pos[2])).rjust(8, " ")

	stupid_number_format = str(i).rjust(7, " ")
	pdb_out.write("ATOM" + stupid_number_format + "  C   GLY     1" + ridiculous_pos_x + ridiculous_pos_y + ridiculous_pos_z + "\n")
	i += 1

pdb_out.write("TER\nENDMDL\n")

pdb_out.close()
