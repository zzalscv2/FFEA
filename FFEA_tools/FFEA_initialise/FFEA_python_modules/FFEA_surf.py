import sys
import numpy as np
import FFEA_triangle

class FFEA_surf:

	def __init__(self, surf_fname):
		
		# Open file and check type
		fin = open(surf_fname, "r")
		if fin.readline().strip() == "surfacemesh":
			fin.close()
			self.load_from_surf(surf_fname)
		else:
			sys.exit("Error. This is not a surface file recognized by FFEA.\n")
	
	def load_from_surf(self, fname):
		
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


	def calc_area(self):
	
		self.area = 0.0
		for f in self.face:
			self.area += f.calc_area(self.node)

		return self.area
