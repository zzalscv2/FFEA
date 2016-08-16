import numpy as np
import sys

class FFEA_vdw:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "VdW file " + fname  + " not found.\nReturning empty object."
			return

		# Header
		if fin.readline().rstrip() != "ffea vdw file":
			print "Error. Expected to read 'ffea vdw file'. This may not be an ffea vdw file"
			return

		# num_faces
		try:
			self.num_faces = int(fin.readline().split()[1])
			self.vdw_index = np.array([-1 for i in range(self.num_faces)])

		except(ValueError):
			print "Error. Expected to read:"
			print "num_faces %d"
			self.reset()
			fin.close()
			return

		# Vdw indices
		if fin.readline().strip() != "vdw params:":
			print "Error. Expected to read 'vdw params:' to begin the vdw indices section."
			self.reset()
			fin.close()
			return

		for i in range(self.num_faces):
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				self.vdw_index[i] = int(line)

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_faces = " + str(self.num_faces) + "\nnum_faces read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected a vdw index of the form '%d' for face " + str(i) + ", but found " + line
				self.reset()
				fin.close()
				return
		
		fin.close()

	def write_to_file(self, fname):

		fout = open(fname, "w")

		# Write header info
		fout.write("ffea vdw file\nnum_faces %d\nvdw params:\n" % (self.num_faces))

		# Write indices
		for ind in self.vdw_index:
			fout.write(str(ind) + "\n")

		fout.close()

	def calc_active_areas(self, surf, node):

		areas = [0.0 for i in range(7)]
		for i in range(self.num_faces):
			areas[self.vdw_index[i] + 1] += surf.face[i].calc_area(node) 
		
		return areas

	def set_num_faces(self, num_faces):
	
		self.num_faces = num_faces
		self.make_inactive()

	def make_inactive(self):

		self.vdw_index = [-1 for i in range(self.num_faces)]
	
	def reset(self):
		self.vdw_index = []
		self.num_faces = 0
