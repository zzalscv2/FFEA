import sys, os
import math

def read_pdb_frame(file):
        while True:
                line = file.readline()
                if line == "":
                        break
                line = line.replace("-", " -")
                if "REMARK" in line:
                        continue
                if "MODEL" in line:
                        continue
                if "ATOM" in line:
                        frame = []
                        while True:
                                if "ENDMDL" in line:
                                        break
                                if "TER" in line:
                                        break
                                if "END" in line:
                                        break
				# ignore the first 31 characters of the line
				line = line[31:]
				sline = line.split()

                                frame.append([float(sline[0]), float(sline[1]), float(sline[2])])
                                line = (file.readline())
                                line = line.replace("-", " -")
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

if len(sys.argv) != 9:
	sys.exit("Usage: python transform_pdb.py [INPUT PDB] [YAW] [PITCH] [ROLL] [DX] [DY] [DZ] [OUTPUT PDB]")

pdb_in = open(sys.argv[1], "r")
num_atoms, frame = read_pdb_frame(pdb_in)

centroid = get_centroid(num_atoms, frame)
print "centroid =", centroid

yaw = float(sys.argv[2])
pitch = float(sys.argv[3])
roll = float(sys.argv[4])

dx = float(sys.argv[5])
dy = float(sys.argv[6])
dz = float(sys.argv[7])

rotation_matrix = construct_rotation_matrix(yaw, pitch, roll)
print "rotation matrix =", rotation_matrix

pdb_out = open(sys.argv[8], "w")

frame = [[x[0] - centroid[0], x[1] - centroid[1], x[2] - centroid[2]] for x in frame]

while True:
	i = 0
	for pos in frame:
		new_pos = mat_mult(rotation_matrix, pos)
		new_pos[0] += centroid[0] + dx
		new_pos[1] += centroid[1] + dy
		new_pos[2] += centroid[2] + dz

		ridiculous_pos_x = ("%.3f" % (new_pos[0])).rjust(12, " ")
		ridiculous_pos_y = ("%.3f" % (new_pos[1])).rjust(8, " ")
		ridiculous_pos_z = ("%.3f" % (new_pos[2])).rjust(8, " ")

		stupid_number_format = str(i).rjust(7, " ")
		pdb_out.write("ATOM" + stupid_number_format + "  C   GLY     1" + ridiculous_pos_x + ridiculous_pos_y + ridiculous_pos_z + "\n")
		i += 1

	pdb_out.write("TER\nENDMDL\n")

	num_atoms, frame = read_pdb_frame(pdb_in)
	if num_atoms == -1:
		break
	frame = [[x[0] - centroid[0], x[1] - centroid[1], x[2] - centroid[2]] for x in frame]

pdb_in.close()
pdb_out.close()
