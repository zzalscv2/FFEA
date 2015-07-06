import sys, os
import numpy as np

class FFEA_kinetic_map:

	def __init__(self, map_fname):
		
		# Check for initialisation type
		if map_fname == "":
			print "No map supplied. Initialising empty map objects..."
			self.dense = None
			self.sparse = None
			self.map_type = "d"
			self.num_rows = 0
			self.num_columns = 0
			self.num_entries = 0
			return

		# Open file and check initial stuff
		fin = open(map_fname, "r")
		line = fin.readline()
		if "(Dense)" in line:
			self.map_type = "d"
		elif "(Sparse)" in line:
			self.map_type = "s"
		else:
			sys.exit("Error. Expected but didn't find 'FFEA Kinetic Conformation Mapping File (Dense/Sparse)'\n")

		self.num_columns = int(fin.readline().split()[1])
		self.num_rows = int(fin.readline().split()[1])
		self.num_entries = int(fin.readline().split()[1])
		fin.readline()
		
		if self.map_type == "d":
			self.dense = np.array([[0.0] * self.num_columns for i in range(self.num_rows)])
			self.sparse = None
			print "Reading dense matrix. No sparse will be created..."
			for i in range(self.num_rows):
				sline = fin.readline().split()
				for j in range(self.num_columns):
					self.dense[i][j] = sline[j]
			print "...done!"

		fin.close()
		
	def set_dense(self, dense_map):

		self.dense = dense_map
		self.map_type = "d"
		self.num_rows = len(self.dense)
		self.num_columns = len(self.dense[0])
		for i in range(self.num_rows):
			for j in range(self.num_columns):
				if abs(self.dense[i][j]) < 0.001:
					self.dense[i][j] = 0
				else:
					self.num_entries += 1

	def write_to_file(self, fname):

		# Open file and write initial stuff
		fout = open(fname, "w")
		fout.write("FFEA Kinetic Conformation Mapping File (Dense)\nnum_nodes_from %d\nnum_nodes_to %d\nnum_entries %d\nmap:\n" % (self.num_columns, self.num_rows, self.num_entries))

		# Write map
		for i in range(self.num_rows):
			for j in range(self.num_columns):
				fout.write("%f " % (self.dense[i][j]))
			fout.write("\n")

		fout.close()
	
	def apply_to_nodes(self, FFEA_nodes):

		# Check if application is possible
		if self.num_columns != FFEA_nodes.num_nodes:
			sys.exit("Error. Matrix requires " + self.num_columns + " nodes. '" + innode + "' only has " + FFEA_nodes.num_nodes + " nodes.\n")

		# Make return object
		out_nodes = np.array([[0.0,0.0,0.0] for i in range(self.num_rows)])

		# Apply to nodes
		for i in range(self.num_rows):
			for j in range(self.num_columns):
				for k in range(3):
					out_nodes[i][k] += self.dense[i][j] * FFEA_nodes.pos[j][k]

		return out_nodes
				
	def apply_to_map(self, FFEAmap):

		# Check if application is possible
		if self.num_columns != FFEAmap.num_rows:
			sys.exit("Error. Matrix cannot be applied.\n")

		# Make return map
		outmap = FFEA_kinetic_map("")
		
		# Apply map
		outmap.dense = np.dot(self.dense, FFEAmap.dense)
		outmap.set_dense(outmap.dense)
		
		return outmap
		
