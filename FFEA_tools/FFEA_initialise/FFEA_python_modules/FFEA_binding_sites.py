import numpy as np
import sys

class FFEA_binding_sites:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. File " + fname  + " not found. Creating empty object...done!"
			self.reset()
			return

		# Header
		if fin.readline().rstrip() != "ffea binding sites file":
			print "Error. Expected to read 'ffea binding sites file'. This may not be an ffea binding sites file"
			return

		# num_binding_sites
		try:
			self.num_binding_sites = int(fin.readline().split()[1])
			self.bsites = [FFEA_binding_site() for i in range(self.num_binding_sites)]

		except(ValueError):
			print "Error. Expected to read:"
			print "num_binding_sites = %d"
			self.reset()
			fin.close()
			return

		# All sites
		if fin.readline().strip() != "binding sites:":
			print "Error. Expected to read 'binding sites:' to begin the binding sites section."
			self.reset()
			fin.close()
			return

		for i in range(self.num_binding_sites):
			indices = []
			num_faces
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				sline = line.split()
				num_faces = int(sline[0])
				if num_faces != len(sline) - 1:
					raise IndexError

				self.bsites[i].set_structure([int(index) for index in sline[1:]])

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_faces = " + str(self.num_faces) + "\nnum_faces read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected " + str(num_faces) + " faces for binding site " + str(i) + ", but found " + str(len(sline) - 1)
				self.reset()
				fin.close()
				return
		
		fin.close()

	def write_to_file(fname)

		fout = open(fname, "w")

		# Write the header info
		fout.write("ffea binding sites file\nnum_binding_sites = %d\nbinding sites:\n" % (self.num_binding_sites))
		
		# Write the sites
		for s in self.bsites:
			fout.write("%d" % (s.num_faces))
			
			for f in s.face_index:
				fout.write(" %d" % (f))
			fout.write("\n")
		fout.close()

	def add_site(self, bsite):

		self.bsites.append(bsite)
		self.num_binding_sites += 1

	def reset(self):
		self.bsites = []
		self.num_binding_sites = 0

class FFEA_binding_site:

	def __init__(self):

		self.reset()

	def reset(self):

		self.num_faces = 0
		self.face_index = []

	def set_structure(self, faces):
		
		self.num_faces = len(faces)
		self.face_index = faces
