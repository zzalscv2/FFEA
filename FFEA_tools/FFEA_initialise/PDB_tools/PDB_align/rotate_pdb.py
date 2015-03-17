import sys, os
import math

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

if len(sys.argv) != 6:
	sys.exit("Usage: python rotate_pdb.py [INPUT PDB] [YAW] [PITCH] [ROLL] [OUTPUT PDB]")

pdb_in = open(sys.argv[1], "r")
num_atoms, frame = read_pdb_frame(pdb_in)
pdb_in.close()

centroid = get_centroid(num_atoms, frame)
print "centroid =", centroid

yaw = math.radians(float(sys.argv[2]))
pitch = math.radians(float(sys.argv[3]))
roll = math.radians(float(sys.argv[4]))

rotation_matrix = construct_rotation_matrix(yaw, pitch, roll)
print "rotation matrix =", rotation_matrix

frame = [[x[0] - centroid[0], x[1] - centroid[1], x[2] - centroid[2]] for x in frame]

pdb_out = open(sys.argv[5], "w")
i = 0
for pos in frame:
	new_pos = mat_mult(rotation_matrix, pos)
#	new_pos[0] += centroid[0]
#	new_pos[1] += centroid[1]
#	new_pos[2] += centroid[2]

	ridiculous_pos_x = ("%.3f" % (new_pos[0])).rjust(12, " ")
	ridiculous_pos_y = ("%.3f" % (new_pos[1])).rjust(8, " ")
	ridiculous_pos_z = ("%.3f" % (new_pos[2])).rjust(8, " ")

	stupid_number_format = str(i).rjust(7, " ")
	pdb_out.write("ATOM" + stupid_number_format + "  C   GLY     1" + ridiculous_pos_x + ridiculous_pos_y + ridiculous_pos_z + "\n")
	i += 1

pdb_out.write("TER\nENDMDL\n")

pdb_out.close()
