import FFEA_node, FFEA_surf, FFEA_vdw, FFEA_pin, FFEA_topology
from Vectors import *
from math import *
import numpy as np

class pyFFEA:
	
	def __init__(self):
		
		num_nodes = 0
		num_interior_nodes = 0
		num_surface_nodes = 0
		nodes = []

		num_surf_faces = 0
		surf = []

		num_elements = 0
		top = []
	
		num_vdw_faces = 0
		vdw = []

		num_pinned_nodes = 0
		pin = []

	def load_script(script_fname):
		
		# Open script
		fin = open(script_fname, "r")

		# Check if script file
		if fin.readline().strip() != "FFEA Script File":
                        sys.exit("Error. '%s' first line is not 'FFEA Script File', so possibly not a script file!\n" % (script_fname))
				
		# First section should be system
		while fin.readline().strip() != "<system>":
			continue

		# Do stuff until end of system!
		line = fin.readline().strip()
		while != "</system>":
			line = line.replace("<", "").replace(">", "")
			lvalue, rvalue = line.split("=")
			lvalue = lvalue.strip()
			rvalue = rvalue.strip()
			if lvalue = "dt":
				self.system_params.dt = float(rvalue)
			if lvalue = "num_steps":
				self.system_params.num_steps = int(rvalue)

	def load_nodes(nodes_fname):

		# Open nodes file
		fin = open(nodes_fname, "r")

		# Check if nodes file
		if fin.readline().strip() != "FFEA Nodes File":
			sys.exit("Error. '%s' first line is not 'FFEA Nodes File', so possibly not a nodes file!\n" % (nodes_fname))

		# Get num_nodes etc
		self.num_nodes = int(fin.readline().split()[1])
		self.num_interior_nodes = int(fin.readline().split()[1])
		self.num_surface_nodes = int(fin.readline().split()[1]) 
 		
		# Get interior nodes
		fin.readline()
		for i in range(self.num_interior_nodes):
			sline = fin.readline().split()
			nodes.append(np.array([float(sline[0]), float(sline[1]), float(sline[2])]), float)

		# Get surface nodes
		fin.readline()
		for i in range(self.num_sirface_nodes):
			sline = fin.readline().split()
			nodes.append(np.array([float(sline[0]), float(sline[1]), float(sline[2])]), float)
		
		# Close file and exit
		fin.close()
