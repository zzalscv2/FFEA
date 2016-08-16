import numpy as np
import sys

class FFEA_stokes:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. Stokes file " + fname  + " not found."
			return

		# Header
		if fin.readline().rstrip() != "ffea stokes radii file":
			print "Error. Expected to read 'ffea stokes radii file'. This may not be an ffea stokes radii file"
			return

		# num_nodes
		try:
			self.num_nodes = int(fin.readline().split()[1])
			self.stokes_radius = np.array([1.0 for i in range(self.num_nodes)])

		except(ValueError):
			print "Error. Expected to read:"
			print "num_nodes = %d"
			self.reset()
			fin.close()
			return

		# stokes indices
		for i in range(self.num_nodes):
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				self.stokes_radius[i] = float(line)

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_nodes = " + str(self.num_nodes) + "\nnum_nodes read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected a stokes radius of the form '%f' for node " + str(i) + ", but found " + line
				self.reset()
				fin.close()
				return
		
		fin.close()

	def reset(self):
		self.stokes_radius = []
		self.num_nodes = 0
