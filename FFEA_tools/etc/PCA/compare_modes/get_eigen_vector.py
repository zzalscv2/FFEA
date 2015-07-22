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

if len(sys.argv) != 3:
	sys.exit("Usage: python get_eigen_vector.py [INPUT PDB] [OUTPUT EIGENVECTOR FILE]")

pdb_in = open(sys.argv[1], "r")
num_atoms, frame1 = read_pdb_frame(pdb_in)
num_atoms, frame2 = read_pdb_frame(pdb_in)
pdb_in.close()

eigen = []
sum = 0
for i in xrange(num_atoms):
	adx = (frame2[i])[0] - (frame1[i])[0]
	ady = (frame2[i])[1] - (frame1[i])[1]
	adz = (frame2[i])[2] - (frame1[i])[2]
	eigen.append(adx)
	eigen.append(ady)
	eigen.append(adz)
	sum += adx * adx + ady * ady + adz * adz
l = math.sqrt(sum)

eigen_out = open(sys.argv[2], "w")
eigen_out.write("N = " + str(len(eigen)) + "\n")
for n in eigen:
	eigen_out.write(str(n/l) + "\n")
eigen_out.close()
