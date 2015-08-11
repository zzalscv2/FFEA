import sys
import numpy as np
import FFEA_triangle

class FFEA_surf:

	def __init__(self, surf_fname):
		
		# Open file and check type
		fin = open(surf_fname, "r")
		line = fin.readline().strip()
		if line == "surfacemesh":
			fin.close()
			self.load_from_netgen_surf(surf_fname)
		elif line == "ffea surface file":
			fin.close()
			self.load_from_ffea_surf(surf_fname)
		else:
			sys.exit("Error. This is not a surface file recognized by FFEA.\n")
	
	def load_from_netgen_surf(self, fname):
		
		fin = open(fname, "r")
		fin.readline()
		
		# Get nodes into numpy array
		self.num_nodes = int(fin.readline().strip())
		temp_nodes = []
		for i in range(self.num_nodes):
			sline = fin.readline().split()
			temp_nodes.append(np.array([float(sline[0].strip()), float(sline[1].strip()), float(sline[2].strip())]))
		
		self.node = np.array(temp_nodes)
		del temp_nodes

		# Now get faces as FFEA_traingles
		self.num_faces = int(fin.readline().strip())
		self.face = []
		for i in range(self.num_faces):
			sline = fin.readline().split()
			self.face.append(FFEA_triangle.FFEA_triangle(int(sline[0].strip()) - 1, int(sline[1].strip()) - 1, int(sline[2].strip()) - 1))

		fin.close()


	def load_from_ffea_surf(self, fname):
		
		fin = open(fname, "r")

		# Ignoring nodes
		self.num_nodes = 0
		self.node = []

		# Get initial stuff
		fin.readline()
		self.num_faces = int(fin.readline().split()[1])
		self.face = []
		fin.readline()
	
		# Now get faces as FFEA_traingles
		for i in range(self.num_faces):
			sline = fin.readline().split()
			self.face.append(FFEA_triangle.FFEA_triangle(int(sline[1].strip()), int(sline[2].strip()), int(sline[3].strip())))

		fin.close()

	def write_to_netgen_surf(self, fname):
	
		fout = open(fname, "w")
		fout.write("surfacemesh\n")
		fout.write(str(self.num_nodes) + "\n")
		for anode in self.node:
			fout.write("%6.3f %6.3f %6.3f\n" % (anode[0], anode[1], anode[2]))

		fout.write(str(self.num_faces) + "\n")
		for aface in self.face:
			fout.write("%d %d %d\n" % (aface.n[0] + 1, aface.n[1] + 1, aface.n[2] + 1))

		fout.close()

	def scale(self, scale):
		for i in range(len(self.node)):
			self.node[i] *= scale

	def calc_area(self):
	
		self.area = 0.0
		for f in self.face:
			self.area += f.calc_area(self.node)

		return self.area

	def calc_face_centroid(self, face_index, node):

		return (node.pos[self.face[face_index].n[0]] + node.pos[self.face[face_index].n[1]] + node.pos[self.face[face_index].n[2]]) / 3.0
		
