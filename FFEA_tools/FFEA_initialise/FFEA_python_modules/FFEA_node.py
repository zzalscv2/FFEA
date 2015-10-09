import numpy as np
import sys

class FFEA_node:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. File " + fname  + " not found."
			return

		# Header
		if fin.readline().rstrip() != "ffea node file":
			print "Error. Expected to read 'ffea node file'. This may not be an ffea node file"
			return

		# num_nodes
		try:
			self.num_nodes = int(fin.readline().split()[1])
			self.num_surface_nodes = int(fin.readline().split()[1])
			self.num_interior_nodes = int(fin.readline().split()[1])

		except(ValueError):
			print "Error. Expected to read:"
			print "num_nodes = %d\nnum_surface_nodes = %d\nnum_interior_nodes = %d"
			self.reset()
			fin.close()
			return

		# Surface nodes
		if fin.readline().strip() != "surface nodes:":
			print "Error. Expected to read 'surface nodes:' to begin the surface node section."
			self.reset()
			fin.close()
			return

		for i in range(self.num_surface_nodes):
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				sline = line.split()
				self.pos.append([float(sline[0]), float(sline[1]), float(sline[2])])

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_nodes = " + str(self.num_nodes) + "\nnum_nodes read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected a node position of the form '%f %f %f' for node " + str(i) + ", but found " + line
				self.reset()
				fin.close()
				return
		
		# Interior nodes
		if fin.readline().strip() != "interior nodes:":
			print "Error. Expected to read 'interior nodes:' to begin the interior node section."
			self.reset()
			fin.close()
			return
		
		for i in range(self.num_surface_nodes, self.num_nodes):
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				sline = line.split()
				self.pos.append([float(sline[0]), float(sline[1]), float(sline[2])])

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_nodes = " + str(self.num_nodes) + "\nnum_nodes read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected a node position of the form '%f %f %f' for node " + str(i) + ", but found " + line
				self.reset()
				fin.close()
				return

		# Convert to numpy
		self.pos = np.array(self.pos)
		fin.close()

	def write_to_file(self, fname):
		
		fout = open(fname, "w")
		fout.write("ffea node file\nnum_nodes %d\nnum_surface_nodes %d\nnum_interior_nodes %d\nsurface nodes:\n" % (self.num_nodes, self.num_surface_nodes, self.num_interior_nodes))
		for i in range(self.num_surface_nodes):
			fout.write("%8.6f %8.6f %8.6f\n" % (self.pos[i][0], self.pos[i][1], self.pos[i][2]))
		fout.write("interior nodes:\n")
		for i in range(self.num_surface_nodes, self.num_nodes):
			fout.write("%8.6f %8.6f %8.6f\n" % (self.pos[i][0], self.pos[i][1], self.pos[i][2]))
		fout.close()

	def reset(self):
		self.pos = []
		self.num_nodes = 0
		self.num_interior_nodes = 0
		self.num_surface_nodes = 0

	def scale(self, scale):
		self.pos *= scale
		
