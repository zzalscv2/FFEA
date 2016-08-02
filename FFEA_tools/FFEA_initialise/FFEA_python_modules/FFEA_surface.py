import os
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
		if not os.path.exists(fname):
			print("\tFile '" + fname + "' not found.")
			return
	
		# File format?
		base, ext = os.path.splitext(fname)
		if ext == ".surf":
			try:
				self.load_surf(fname)
				self.valid = True
			except:
				print("\tUnable to load FFEA_surface from " + fname + ". Returning empty object...")

		elif ext == ".face":
			try:
				self.load_face(fname)
				self.valid = True
			except:
				print("\tUnable to load FFEA_surface from " + fname + ". Returning empty object...")

		elif ext == ".vol":
			try:
				self.load_vol(fname)
				self.valid = True
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

	def load_face(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		sline = fin.readline().split()
		if len(sline) != 2:
			print("\tExpected '<num_faces> 1' but found " + line)
			raise TypeError

		num_faces = int(sline[0])

		# Read faces now	
		while(True):
			sline = fin.readline().split()

			if sline[0].strip() == "#":
				break

			# Get a face
			sline = sline[1:4]
			f = FFEA_face_tri_lin()
			f.set_indices(sline)
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
	
	def write_to_file(self, fname):

		print "Writing to " + fname + "..."

		# Write differently depending on format
		base, ext = os.path.splitext(fname)

		if ext == ".vol":
			fout = open(fname, "a")
			fout.write("# surfnr    bcnr   domin  domout      np      p1      p2      p3\nsurfaceelementsgi\n%d\n" % (self.num_faces))
			findex = 1
			for f in self.face:
				#findex += 1
				fout.write(" %d 1 1 0 %d" % (findex, len(f.n)))
				for n in f.n:
					fout.write(" %d" % (n))

				fout.write("\n")

			fout.write("\n\n")
		fout.close()
		print "done!"

	def reset(self):

		self.face = []
		self.num_faces = 0
		self.valid = False

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

	def isSame(self, face):
		n = np.array(self.n)
		m = np.array(face.n[0:3])

		# Do all permutations of indices
		cycle = np.array([[0,0,1],[1,0,0],[0,1,0]])
		swap = np.array([[0,1,0],[1,0,0],[0,0,1]])

		# The faces should be in the opposite order (for connectivity reasons), so swap first
		m = np.dot(swap, m)

		for i in range(3):
			m = np.dot(cycle, m)
			if np.array_equal(n, m):
				return True

		# If not in this order, then mesh is wrong anyway, so return False
		#m = np.dot(swap, m)
		#for i in range(3):
		#	m = np.dot(cycle, m)
		#	print n, m
		#	if np.array_equal(n, m):
		#		print "Equal!"
		#		return True
		
		return False

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
