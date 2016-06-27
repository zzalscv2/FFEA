from os import path
from time import sleep
import numpy as np

class FFEA_surface:

	def __init__(self, fname = ""):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname):

		print("Loading FFEA surface file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".surf":
			try:
				self.load_surf(fname)
			except:
				print("\tUnable to load FFEA_surface from " + fname + ". Returning empty object...")

		elif ext == ".vol":
			try:
				self.load_vol(fname)
			except:
				print("\tUnable to load FFEA_surface from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_surf(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea surface file" and line != "walrus surface file":
			print("\tExpected 'ffea surf file' but found " + line)
			raise TypeError

		num_faces = int(fin.readline().split()[1])

		fin.readline()

		# Read faces now	
		while(True):
			sline = fin.readline().split()

			if len(sline) == 0:
				break

			# Get a face	
			if len(sline) == 3 or len(sline) == 4:
				f = FFEA_face_tri_lin()

				if len(sline) == 3:
					# Just a face
					f.set_indices(sline)
				else:
					# Face with parent element
					f.set_indices(sline[1:], elindex = sline[0])

			elif len(sline) == 6 or len(sline) == 7:
				f = FFEA_face_tri_sec()

				if len(sline) == 6:
					# Just a face
					f.set_indices(sline)
				else:
					# Face with parent element
					f.set_indices(sline[1:], elindex = sline[0])

			self.add_face(f)

		fin.close()

	def add_face(self, f):

		self.face.append(f)
		self.num_faces += 1
		
	def print_details(self):

		print "num_faces = %d" % (self.num_faces)
		sleep(1)

		index = -1
		for f in self.face:
			index += 1
			outline = "Face " + str(index) + ": "

			for n in f.n:
				outline += str(n) + " "
			
			if f.elindex != -1:
				outline += ", Parent Element = %d" % (f.elindex)

			print outline
	
	def reset(self):

		self.face = []
		self.num_faces = 0

class FFEA_face:

	def __init__(self):

		self.reset()

	def set_indices(self, alist, elindex = -1):

		# Test for correct number of nodes
		if len(alist) != len(self.n):
			print "Incorrect number of nodes for assignment to this face type."
			return

		for i in range(len(self.n)):
			self.n[i] = int(alist[i])

		self.elindex = int(elindex)

	def calc_centroid(self, node):
	
		centroid = np.array([0.0,0.0,0.0])
		for i in self.n:
			centroid += node.pos[i]
			
		return centroid * (1.0 / len(self.n))

	def reset(self):
		
		self.n = []
		self.elindex = -1

class FFEA_face_tri_lin(FFEA_face):

	def reset(self):

		self.n = [0,1,2]
		self.elindex = -1

class FFEA_face_tri_sec(FFEA_face):

	def reset(self):

		self.n = [0,1,2,3,4,5]
		self.elindex = -1
