import os, sys
from math import *

class Node:
	def __init__(self, x, y, z, index):
		self.x = x
		self.y = y
		self.z = z
		self.index = index

def process(n, rmsd, rmsd_nodes, vector_x, vector_y, vector_z, centroid_x, centroid_y, centroid_z):
	line = inputtraj.readline()
	if line == "STATIC\n":
		return 0.0
	else:
		num_nodes = int(inputtraj.readline())
		dr_x = 0.0
		dr_y = 0.0
		dr_z = 0.0
		for i in range(len(rmsd_nodes)):
			if i == 0:
				for j in range(rmsd_nodes[i].index + 1):
					line = inputtraj.readline()
			else:
				for j in range(rmsd_nodes[i].index - rmsd_nodes[i-1].index):
					line = inputtraj.readline()
			if n == 0:
				sline = line.split()
				dr_x = dr_x + float(sline[0]) * 1e10 
				dr_y = dr_y + float(sline[0]) * 1e10
				dr_z = dr_z + float(sline[0]) * 1e10 
		
		#Finishing Reading Blob
		for i in range(num_nodes - (rmsd_nodes[-1].index + 1)):
			line = inputtraj.readline()

		dr_x = dr_x / len(rmsd_nodes) - centroid_x
		dr_y = dr_y / len(rmsd_nodes) - centroid_y
		dr_z = dr_z / len(rmsd_nodes) - centroid_z

		#dot product
		dot = dr_x * vector_x + dr_y * vector_y + dr_z * vector_z	
	return pow(dot, 2)
			 
	
if len(sys.argv) != 8:
	sys.exit("Usage: python rmsd_radially.py [INPUT .NODE FILE] [INPUT TRAJ FILE] [BLOB NUMBER] [TOTAL NUMBER OF BLOBS] [NUMBER OF FRAMES] [LOWER LIMIT RADIUS (angstroms)] [UPPER LIMIT RADIUS (angstroms)]")

blob_number = sys.argv[3]
num_blobs = int(sys.argv[4])
num_frames = int(sys.argv[5])
lower_radius = float(sys.argv[6])
upper_radius = float(sys.argv[7])

num_nodes = [0,0,0]

#Finding relevant nodes
inputnode = open(sys.argv[1], "r")
inputnode.readline() #header
num_nodes = int(inputnode.readline().split()[1]) #num_nodes
num_surface_nodes = int(inputnode.readline().split()[1]) #num_surface_nodes
num_interior_nodes = int(inputnode.readline().split()[1]) #num_interior_nodes
inputnode.readline() #"surface nodes:"

rmsd_nodes = []

for i in range(num_surface_nodes):
	sline = inputnode.readline().split()
	x_pos = float(sline[0]) * 1e-10
	y_pos = float(sline[1]) * 1e-10
	z_pos = float(sline[2]) * 1e-10
	radius = sqrt(pow(x_pos, 2) + pow(y_pos, 2) + pow(z_pos, 2))
	if radius > lower_radius and radius < upper_radius:
		 rmsd_nodes.append(Node(x_pos, y_pos, z_pos, i))

inputnode.readline() #"interior nodes:"
for i in range(num_interior_nodes):
	sline = inputnode.readline().split()
	x_pos = float(sline[0])
	y_pos = float(sline[1])
	z_pos = float(sline[2])
	radius = sqrt(pow(x_pos, 2) + pow(y_pos, 2) + pow(z_pos, 2))
	if radius > lower_radius and radius < upper_radius:
		  rmsd_nodes.append(Node(x_pos, y_pos, z_pos, i + num_surface_nodes))
inputnode.close()
print "Calculating <(x-x_0)^2> for " + str(len(rmsd_nodes)) + " out of " + str(num_nodes) + " nodes\n"

#Calculating and normalising directional vector
vector_x = 0.0
vector_y = 0.0
vector_z = 0.0
for node in rmsd_nodes:
	vector_x = vector_x + node.x
	vector_y = vector_y + node.y
	vector_z = vector_z + node.z

vector_x = vector_x / len(rmsd_nodes)
vector_y = vector_y / len(rmsd_nodes)
vector_z = vector_z / len(rmsd_nodes)
centroid_x = vector_x
centroid_y = vector_y
centroid_z = vector_z
mod = sqrt(pow(vector_x, 2) + pow(vector_y, 2) + pow(vector_z, 2))
vector_x = vector_x / mod
vector_y = vector_y / mod
vector_z = vector_z / mod
print "Centroid of selected nodes: (" + str(centroid_x) + ", " + str(centroid_y) + ", " + str(centroid_z) + ")\n"

#Calculating rmsd of section of blob
inputtraj = open(sys.argv[2], "r")
rmsd = 0.0

for i in range(num_frames):
	inputtraj.readline() #*
	for blob in range(num_blobs):
		sline = inputtraj.readline().split()
		
		if sline[1][0] == blob_number:
			if sline[3] == 0 or sline[3] == 1:
				rmsd = rmsd + process(1, rmsd, rmsd_nodes, vector_x, vector_y, vector_z, centroid_x, centroid_y, centroid_z)	
			else:
				rmsd = rmsd + process(0, rmsd, rmsd_nodes, vector_x, vector_y, vector_z, centroid_x, centroid_y, centroid_z)
		else:
			rmsd = rmsd + process(1, rmsd, rmsd_nodes, vector_x, vector_y, vector_z, centroid_x, centroid_y, centroid_z)
	
rmsd = rmsd / float(num_frames) 
print "<(x-x_0)^2> = " + str(rmsd) + "Angstroms^2 for your selected section and frame length."	
