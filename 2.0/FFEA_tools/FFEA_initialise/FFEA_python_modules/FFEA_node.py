import sys, os
import numpy as np

class FFEA_node:

	def __init__(node_fname):
		
		# Open file and check initial stuff
		fin = open(node_fname, "r")
		if(fin.readline().strip() != "ffea node file"):
			sys.exit("Error. Expected but didn't find 'ffea node fname'\n")
		
		self.num_nodes = int(fin.readline().split()[1].strip())
		self.num_surface_nodes = int(fin.readline().split()[1].strip())
		self.num_interior_nodes = int(fin.readline().split()[1].strip())
		self.pos = np.array([0.0, 0.0, 0.0] for i in range(self.num_nodes))

		# Surface nodes
		fin.readline()
		for i in range(self.num_surface_nodes):
			sline = fin.readline().split()
			self.pos[i][0] = float(sline[0])
			self.pos[i][1] = float(sline[1])
			self.pos[i][2] = float(sline[2])

		# Interior nodes
                fin.readline()
                for i in range(self.num_surface_nodes, self.num_nodes, 1):
                        sline = fin.readline().split()
                        self.pos[i][0] = float(sline[0])
                        self.pos[i][1] = float(sline[1])
                        self.pos[i][2] = float(sline[2])

		fin.close()
