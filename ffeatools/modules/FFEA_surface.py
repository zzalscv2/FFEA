import os
from time import sleep
import numpy as np

class FFEA_surface:

	def __init__(self, fname = ""):
	
		self.reset()

		self.load(fname)

	def load(self, fname=""):

		print("Loading FFEA surface file...")

		if fname=="":
			print("\tCreating empty object...")
			return
	
		# File format?
		base, ext = os.path.splitext(fname)
		if ext == ".surf":
			self.load_surf(fname)
		elif ext == ".face":
			self.load_face(fname)
		elif ext == ".vol":
			self.load_vol(fname)
		else:
			raise IOError("Unrecognised file extension '" + ext + "'.")

		self.valid = True
		return

	def load_surf(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise IOError("File '" + fname + "' not found.")

		# Test format
		line = fin.readline().strip()
		if line != "ffea surface file" and line != "walrus surface file":
			raise TypeError("Expected 'ffea surf file' but found " + line)

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
			sline = [int(s) - 1 for s in sline]
			f = FFEA_face_tri_lin()
			f.set_indices(sline)
			self.add_face(f)

		fin.close()

	def load_vol(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		lines = fin.readlines()
		fin.close()

		# Test format and get to write place
		i = 0
		line = lines[i].strip()
		while line != "surfaceelements" and line != "surfaceelementsgi":
			i += 1
			try:
				line = lines[i].strip()
			except(IndexError):
				print("\tCouldn't find 'surfaceelements' line. File '" + fname + "' not formatted correctly.")
				self.reset()
				return

			continue

		# Get num_faces
		i += 1
		num_faces = int(lines[i])
		for j in range(i + 1, i + 1 + num_faces):
			try:
				sline = lines[j].split()
				face = FFEA_face_tri_lin()
				face.set_indices([int(sline[5]) - 1, int(sline[6]) - 1, int(sline[7]) - 1])
				self.add_face(face)

			except:
				print("\tCouldn't find the specified %d faces. Only found %d. File '" + fname + "' not formatted correctly." % (num_faces, i))
				self.reset()
				return

	def add_face(self, f):
		self.face.append(f)
		self.num_faces += 1

	def get_element_indices(self, top):

		elcheck = range(top.num_elements)
		index = -1
		for f in self.face:
			index += 1			
			success = 0
			for i in elcheck:
				print "Element ", i, " = ", top.element[i].n
				compare = 0
				for j in range(3):
					if f.n[j] in top.element[i].n:
						compare += 1
				if compare == 3:
					success = 1
					print "Face found = ", f.n
					print "Element removed = ", top.element[i].n
					f.elindex = i
					#elcheck.remove(i)	# Can't remove, as elements can have more than one surface face (especially in coarse structures :( )
					break

			if success != 1:
				print("Face " + str(index) + " could not be found in the topology structure. This topology cannot be paired with this surface. Regenerating surface...")
				print "Face in error = ", self.face[index].n
				for f in self.face:
					f.elindex = None
				return -1

	def upgrade_face(self, index):

		# Replace element with a higher order one
		f = FFEA_face_tri_sec()
		f.n = self.face[index].n
		f.elindex = self.face[index].elindex

		self.face.insert(index, f)
		self.face.pop(index + 1)

	def split_face(self, index):

		# Split second order face into 4 1st orders
		if isinstance(self.face[index], FFEA_face_tri_sec):
			
			# Order is lin0, lin1, lin2, sec01, sec02, sec12
			# So, 4 tris are 034, 153, 245, 354
			oldf = self.face[index]
			f = [FFEA_face_tri_lin() for i in range(4)]
			f[0].set_indices([oldf.n[0], oldf.n[3], oldf.n[4]], oldf.elindex)
			f[1].set_indices([oldf.n[1], oldf.n[5], oldf.n[3]], oldf.elindex)
			f[2].set_indices([oldf.n[2], oldf.n[4], oldf.n[5]], oldf.elindex)
			f[3].set_indices([oldf.n[3], oldf.n[5], oldf.n[4]], oldf.elindex)
			
			# Insert in this order after the old face
			for face in f:
				self.face.append(face)
			
			self.face.pop(index)

			# 1 -> 4 means += 3
			self.num_faces += 3

	def check_normals(self, node, top):
		
		# Get element for each face and make normal point away from it
		for i in range(self.num_faces):
			elindex = self.face[i].elindex
			
			# El to face vector
			cf = self.face[i].calc_centroid(node) - top.element[elindex].calc_centroid(node)

			# Current face normal
			norm = np.cross(node.pos[self.face[i].n[1]] - node.pos[self.face[i].n[0]], node.pos[self.face[i].n[2]] - node.pos[self.face[i].n[0]])

			# Switch if in different directions
			if np.dot(cf, norm) < 0.0:
				index = [self.face[i].n[1], self.face[i].n[2]]
				self.face[i].n[1] = index[1]
				self.face[i].n[2] = index[0]
				for j in range(4):
					if top.element[elindex].n[j] == index[0]:
						top.element[elindex].n[j] = index[1]
					elif top.element[elindex].n[j] == index[1]:
						top.element[elindex].n[j] = index[0]

	def print_details(self):

		print "num_faces = %d" % (self.num_faces)
		sleep(1)

		index = -1
		for f in self.face:
			index += 1
			outline = "Face " + str(index) + ": "

			for n in f.n:
				outline += str(n) + " "
			
			if f.elindex != None:
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
					fout.write(" %d" % (n + 1))

				fout.write("\n")

			fout.write("\n\n")

		elif ext == ".surf":
			fout = open(fname, "w")
			fout.write("ffea surface file\nnum_surface_faces %d\n" % (self.num_faces))
			fout.write("faces:\n")
			for i in range(self.num_faces):
				fout.write("%d " % (self.face[i].elindex))
				for n in self.face[i].n:
					fout.write("%d " % (n))
				fout.write("\n")
		fout.close()
		print "done!"

	def reset(self):

		self.face = []
		self.num_faces = 0
		self.valid = False

class FFEA_face:

	def __init__(self):

		self.reset()

	def set_indices(self, alist, elindex = None):

		# Test for correct number of nodes
		if len(alist) != len(self.n):
			print "Incorrect number of nodes for assignment to this face type."
			return

		for i in range(len(self.n)):
			self.n[i] = int(alist[i])

		try:
			self.elindex = int(elindex)
		except:
			self.elindex = None

	def calc_centroid(self, node):
	
		centroid = np.array([0.0,0.0,0.0])
		for i in self.n:
			centroid += node.pos[i]
			
		return centroid * (1.0 / len(self.n))

	def get_centroid(self):
		return self.centroid

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
		self.elindex = None

class FFEA_face_tri_lin(FFEA_face):

	def reset(self):

		self.n = [0,1,2]
		self.elindex = None

class FFEA_face_tri_sec(FFEA_face):

	def reset(self):

		self.n = [0,1,2,3,4,5]
		self.elindex = None
