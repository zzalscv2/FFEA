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
			print "File " + fname  + " not found.\nCreating empty object...done!"
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
			print "num_binding_sites %d"
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
			try:

				# Shouldn't be at the end yet!
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				# First line (type and stuff)
				sline = line.split()
				site_type = int(sline[1])
				num_faces = int(sline[3])
				
				# Second line is face list
				line = fin.readline()

				# Shouldn't be at the end yet!
				if line == [] or line == None or line == "":
					raise EOFError

				sline = line.split()[1:]
				if(len(sline) != num_faces):
					raise IndexError

				self.bsites[i].set_type(site_type)
				self.bsites[i].set_structure([int(index) for index in sline])

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_faces = " + str(self.num_faces) + "\nnum_faces read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected " + str(num_faces) + " faces for binding site " + str(i) + ", but found " + str(len(sline) - 2)
				self.reset()
				fin.close()
				return
		
		fin.close()

	def write_to_file(self, fname):

		fout = open(fname, "w")

		# Write the header info
		fout.write("ffea binding sites file\nnum_binding_sites %d\nbinding sites:\n" % (self.num_binding_sites))
		
		# Write the sites
		for s in self.bsites:
			fout.write("type %d num_faces %d\nfaces:" % (s.site_type, s.num_faces))
			
			for f in s.face_index:
				fout.write(" %d" % (f))
			fout.write("\n")
		fout.close()

	def add_site(self, bsite):

		if bsite.num_faces == 0:
			print("Your binding site has no faces i.e. is empty. Will not add to structure.\n")
		else:
			self.bsites.append(bsite)
			self.num_binding_sites += 1

	def reset(self):
		self.bsites = []
		self.num_binding_sites = 0

class FFEA_binding_site:

	def __init__(self):

		self.reset()

	def reset(self):

		self.site_type = -1
		self.num_faces = 0
		self.face_index = []

	def set_type(self, site_type):

		self.site_type = site_type

	def set_structure(self, faces):
		
		self.num_faces = len(faces)
		self.face_index = faces
