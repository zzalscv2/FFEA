import sys, os
import numpy as np

class FFEA_topology:
	
	def __init__(top_fname):

		# Open file and check stuff
		fin = open(top_fname, "r")
		if(fin.readline().strip() != "ffea topology file"):
			sys.exit("Error. Expected but didn't find 'ffea topology file.\n")

		self.num_elements = int(fin.readline().split()[1])
		self.num_surface_elements = int(fin.readline().split()[1])
		self.num_interior_elements = int(fin.readline().split()[1])
		self.element = np.array([FFEA_tetra_element(0,0,0,0)] for i in range(num_elements))
		
		# Surface elements
		fin.readline()
		for i in range(self.num_surface_elements):
			sline = fin.readline().split()
			self.element[i].set_nodes(int(sline[0]), int(sline[1]), int(sline[2]), int(sline[3]))	

		# Interior elements
                fin.readline()
                for i in range(self.num_surface_elements, self.num_elements, i):
                        sline = fin.readline().split()
                        self.element[i].set_nodes(int(sline[0]), int(sline[1]), int(sline[2]), int(sline[3]))

		fin.close()

		def calc_volume(FFEA_nodes):

			volume = 0.0
			for element in self.elements:
				volume += element.calc_volume(FFEA_nodes)

			return volume

class FFEA_element:

	def __init__(a,b,c,d):
		
		self.node = np.array([a,b,c,d])
		self.volume = 0.0
		self.surface_area = 0.0
		
	def set_nodes(a,b,c,d):
		
		self.node[0] = a
		self.node[1] = b
		self.node[2] = c
		self.node[3] = d

	def calc_volume(FFEA_nodes):
	
		v = np.array([FFEA_nodes[self.node[i]] - FFEA_nodes[self.node[0]]] for i in range(1,4,1))
		self.volume = np.absolute(np.dot(v[0], np.cross(v[1], v[2]))) / 6.0
		
		return volume
