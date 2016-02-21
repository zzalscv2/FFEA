import numpy as np
import sys

class FFEA_lj:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. LJ file " + fname  + " not found."
			return

		# Header
		if fin.readline().rstrip() != "ffea vdw forcefield params file":
			print "Error. Expected to read 'ffea vdw forcefield params file'. This may not be an ffea vdw forcefield params file"
			return

		# num_faces
		try:
			self.num_face_types = int(fin.readline().split()[1])
			self.lj = [[None for j in range(self.num_face_types)] for i in range(self.num_face_types)]

		except(ValueError):
			print "Error. Expected to read:"
			print "num_vdw_face_types %d"
			self.reset()
			fin.close()
			return

		# lj params
		for i in range(self.num_face_types):
			line = fin.readline()
			if line == "" or line == None:
				print "Error. Expected %d rows but found %d" % (self.num_face_types, i)
				self.reset()
				return
			
			# Get tags
			sline = line.strip().split(") (")
			sline[0] = sline[0][1:]
			sline[-1] = sline[-1][:-1]
			for j in range(self.num_face_types):
				try:
					vdw_eps, r = sline[j].strip().split(",")
					self.lj[i][j] = FFEA_lj_pair(float(vdw_eps), float(r))

				except(IndexError):
					print "Error. Expected %d columns but found %d" % (self.num_face_types, j)
					self.reset()
					return

				except(ValueError):
					print "Error. Expected a block of the form '(%f,%f)' but found " + sline[j]
					self.reset()
					return


	def reset(self):

		self.num_face_types = 0
		self.lj = []

class FFEA_lj_pair:

	def __init__(self, eps, r):
		
		self.epsilon = eps
		self.equilibrium_r = r
