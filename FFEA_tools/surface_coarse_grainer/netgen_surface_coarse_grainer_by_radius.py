#!/usr/bin/env python

import os, sys
import math

class vector3:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	def mag(self):
		return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

class Node:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z


class Face:
	def __init__(self):
		self.smallest_length = 0.0
		self.smallest_length_nodes = [0, 0]
		self.face_nodes = [0,0,0]

	def __init__(self, n0, n1, n2):
		self.smallest_length = 0.0
		self.smallest_length_nodes = [0, 0]
		self.face_nodes = [n0, n1, n2]

	def calc_smallest_length(self, nodes):
		vec01 = vector3(nodes[self.face_nodes[1]].x - nodes[self.face_nodes[0]].x, nodes[self.face_nodes[1]].y - nodes[self.face_nodes[0]].y, nodes[self.face_nodes[1]].z - nodes[self.face_nodes[0]].z) 
		vec02 = vector3(nodes[self.face_nodes[2]].x - nodes[self.face_nodes[0]].x, nodes[self.face_nodes[2]].y - nodes[self.face_nodes[0]].y, nodes[self.face_nodes[2]].z - nodes[self.face_nodes[0]].z) 
		vec12 = vector3(nodes[self.face_nodes[2]].x - nodes[self.face_nodes[1]].x, nodes[self.face_nodes[2]].y - nodes[self.face_nodes[1]].y, nodes[self.face_nodes[2]].z - nodes[self.face_nodes[1]].z) 
		length01 = vec01.mag()
		length02 = vec02.mag()
		length12 = vec12.mag()
		if length01 < length02 and length01 < length12:
			self.smallest_length = length01
			self.smallest_length_nodes = [self.face_nodes[0], self.face_nodes[1]]
		elif length02 < length01 and length02 < length12:
			self.smallest_length = length02
			self.smallest_length_nodes = [self.face_nodes[0], self.face_nodes[2]]
		else:
			self.smallest_length = length12
			self.smallest_length_nodes = [self.face_nodes[1], self.face_nodes[2]]

		
#global functions
# get nodes and data info from input file  
def extract_surface_from_file(infile, num_nodes, num_faces, nodes, faces):
	num_nodes = int(infile.readline())
	for x in range(0, num_nodes):
		line = infile.readline().split()
		nodes.append(Node(float(line[0]), float(line[1]), float(line[2])))
			
	num_faces = int(infile.readline())
	for x in range(0, num_faces):
		line = infile.readline().split()
		faces.append(Face(int(line[0]) - 1, int(line[1]) - 1, int(line[2]) - 1))

# perform coarsening			
def coarsen_surface(nodes, faces, new_nodes, new_faces, threshold, min_radius, max_radius):
	runs = 0
	while True:
		num_faces = len(faces)
		num_nodes = len(nodes)
		print "Run ", runs, "  Nodes ", num_nodes, "  Faces ", num_faces
		global_smallest_length = threshold	
		nodes_involved = [-1, -1]
		for face in faces:
			face.calc_smallest_length(nodes)
			if face.smallest_length < global_smallest_length:
				# See if nodes are within range
				check = 0
				for node in face.smallest_length_nodes:
					this_radius = math.sqrt(math.pow(nodes[node].x, 2) + math.pow(nodes[node].z, 2) + math.pow(nodes[node].z, 2))
					if this_radius > min_radius and this_radius < max_radius:
						check = 1
						break

				if check == 1:
					global_smallest_length = face.smallest_length
					nodes_involved = face.smallest_length_nodes
				
				
		print "Smallest Length ", global_smallest_length, "  Nodes Involved ", nodes_involved[0], nodes_involved[1]	
		if nodes_involved[0] == -1:
			if runs == 0:
				sys.exit("System already coarsened to this length. Pick a larger length or be happy with what God gave you.\n")
			else:
				print "System Coarsened. Ready for output"
				return
	
		# creating a new node at the half way point of this smallest_distance
		new_node = Node((nodes[nodes_involved[0]].x + nodes[nodes_involved[1]].x)/2.0, (nodes[nodes_involved[0]].y + nodes[nodes_involved[1]].y)/2.0, (nodes[nodes_involved[0]].z + nodes[nodes_involved[1]].z)/2.0)
		#print nodes[nodes_involved[0]].x, nodes[nodes_involved[0]].y, nodes[nodes_involved[0]].z
		#print nodes[nodes_involved[1]].x, nodes[nodes_involved[1]].y, nodes[nodes_involved[1]].z
		#print new_node.x, new_node.y, new_node.z
		
		# replacing old line with new node at midpoint
		if nodes_involved[0] > nodes_involved[1]:
			nodes.pop(nodes_involved[0])
			nodes.pop(nodes_involved[1])
		else:
			nodes.pop(nodes_involved[1])
			nodes.pop(nodes_involved[0])

		nodes.append(new_node)
		delete = []
		for face in faces:
			count = 0
			for face_node in face.face_nodes:
				if face_node == nodes_involved[0] or face_node == nodes_involved[1]:
					face.face_nodes[face.face_nodes.index(face_node)] = -1
					count = count + 1
					continue
			#print face.face_nodes
			if count == 2:
				delete.append(faces.index(face))
				#print "Deleting", delete

		for i in reversed(delete):
			faces.pop(i)

		for face in faces:
			change1 = []
			change2 = []
			for face_node in face.face_nodes:	
				if face_node > nodes_involved[0] and face_node > nodes_involved[1]:
					change2.append(face.face_nodes.index(face_node))
					continue
				if face_node > nodes_involved[0] or face_node > nodes_involved[1]:
					change1.append(face.face_nodes.index(face_node))
					continue
		
			for i in change2:
				face.face_nodes[i] = face.face_nodes[i] - 2

			for i in change1:
				face.face_nodes[i] = face.face_nodes[i] - 1
			
			#print face.face_nodes
		for face in faces:
			for face_node in face.face_nodes:
				if face.face_nodes[face.face_nodes.index(face_node)] == -1:
					face.face_nodes[face.face_nodes.index(face_node)] = len(nodes) - 1
			#print face.face_nodes
		runs = runs + 1		
		
					
		
	
if len(sys.argv) != 6:
	sys.exit("Usage: python netgen_surface_coarse_grainer [INPUT .SURF FILE] [OUTPUT .SURF FILE] [LOWER LENGTH THRESHOLD (ANGSTROMS)] [MIN RADIUS (ANGSTROMS)] [MAX RADIUS (ANGSTROMS)]")

infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")
length_thresh = float(sys.argv[3])
min_radius = float(sys.argv[4])
max_radius = float(sys.argv[5])

#testing input file
testline = infile.readline()
if testline[0:11] != "surfacemesh":
	sys.exit("Error. Incorrect format. This isn't a Netgen .surf file\n")

#extracting data from input file
num_nodes = 0
num_faces = 0
nodes = []
faces = []
print "Extracting surface from file"
extract_surface_from_file(infile, num_nodes, num_faces, nodes, faces)
new_num_nodes = num_nodes
new_num_faces = num_faces
new_nodes = []
new_faces = []
print "Coarsening surface to desired level..."
coarsen_surface(nodes, faces, new_nodes, new_faces, length_thresh, min_radius, max_radius)
print "Outputting as Netgen .surf file"
num_nodes = len(nodes)
num_faces = len(faces)
outfile.write("surfacemesh\n")
outfile.write(str(num_nodes) + "\n")
for node in nodes:
	outfile.write(str(node.x) + " " + str(node.y) + " " + str(node.z) + "\n")

outfile.write(str(num_faces) + "\n")
for face in faces:
	outfile.write("\t" + str(face.face_nodes[0] + 1) + "\t" + str(face.face_nodes[1] + 1) + "\t" + str(face.face_nodes[2] + 1) + "\n")

	
