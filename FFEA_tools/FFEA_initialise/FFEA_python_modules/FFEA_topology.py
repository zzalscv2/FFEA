import os
from time import sleep
import numpy as np
import FFEA_surface

class FFEA_topology:

	def __init__(self, fname = ""):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname):

		print("Loading FFEA topology file...")

		# Test file exists
		if not os.path.exists(fname):
			print("\tFile '" + fname + "' not found.")
			return
	
		# File format?
		base, ext = os.path.splitext(fname)
		if ext == ".top":
			try:
				self.load_top(fname)
				self.valid = True
			except:
				print("\tUnable to load FFEA_topology from " + fname + ". Returning empty object...")

		elif ext == ".ele":
			try:
				self.load_ele(fname)
				self.valid = True
			except:
				print("\tUnable to load FFEA_topology from " + fname + ". Returning empty object...")

		elif ext == ".vol":
			try:
				self.load_vol(fname)
				self.valid = True
			except:
				print("\tUnable to load FFEA_topology from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_top(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea topology file" and line != "walrus topology file":
			print("\tExpected 'ffea topology file' but found " + line)
			raise TypeError

		num_elements = int(fin.readline().split()[1])
		num_surface_elements = int(fin.readline().split()[1])
		num_interior_elements = int(fin.readline().split()[1])

		fin.readline()

		# Read elements now
		eltype = 0	
		while(True):
			sline = fin.readline().split()

			# Get an element
			try:
				if sline[0].strip() == "interior":
					eltype = 1
					continue
	
				elif len(sline) == 4:
					el = FFEA_element_tet_lin()
				elif len(sline) == 10:
					el = FFEA_element_tet_sec()

			except(IndexError):
				break

			el.set_indices(sline)
			self.add_element(el)

		fin.close()

	def load_vol(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Find start of elements
		while(True):
			line = fin.readline()
			if line == "":
				raise IOError
			
			elif line.strip() == "volumeelements":
				break

		# Read num_elements
		num_elements = int(fin.readline())

		# Get all elements
		for i in range(num_elements):
			sline = fin.readline().split()[2:]
			el = FFEA_element_tet_lin()
			el.set_indices(sline)
			self.add_element(el)

		# Indexing from 0
		for i in range(self.num_elements):
			for j in range(4):
				self.element[i].n[j] -= 1

	def load_ele(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		sline = fin.readline().split()
		if len(sline) != 3:
			print("\tExpected '<num_elements> <num_vertices> 0' but found " + line)
			raise TypeError

		num_elements = int(sline[0])
		num_interior_elements = 0
		num_surface_elements = num_elements

		# Read elements now
		while(True):
			sline = fin.readline().split()
			if sline[0].strip() == "#":
				break

			sline = sline[1:]
			# Get an element
			try:
				if len(sline) == 4:
					el = FFEA_element_tet_lin()
				elif len(sline) == 10:
					el = FFEA_element_tet_sec()

			except(IndexError):
				break

			el.set_indices(sline)
			self.add_element(el)

		fin.close()

	def add_element(self, el, eltype = -1):

		if eltype == -1:
			self.num_surface_elements += 1
			el.interior = None

		elif eltype == 0:
			self.num_surface_elements += 1
			el.interior = False
		else:
			self.num_interior_elements += 1
			el.interior = True

		self.element.append(el)
		self.num_elements += 1

	def get_num_elements(self):

		return len(self.element)

	def get_linear_nodes(self):
	
		# Get them all
		n = []
		for e in self.element:
			for index in e.n[0:4]:
				n.append(index)
		
		# Make a list of a set
		return list(set(n))
	
	def isElementInterior(self, index):
		
		try:
			testEl = self.element[index]
		except(IndexError):
			print "Element ", index, " does not exist."
			return False

		# First, see if already calculated. Else, set default assumption
		if testEl.interior != None:
			return testEl.interior

		else:
			testEl.interior = False

		# To test if interior, see if a faces are repeated. If all are, element is interior
		i = -1
		faces_not_found = [0,1,2,3]
		for el in self.element:
			i += 1
			el_connected = False
			faces_to_remove = -1

			# Same element doesn't count
			if i == index:
				continue
		
			for j in faces_not_found:
				testFace = testEl.get_linear_face(j)
				
				for k in range(4):
					face = el.get_linear_face(k)

					if face.isSame(testFace):
						faces_to_remove = j
						print j
						el_connected = True
						break

				if el_connected:
					break
	
			try:
				faces_not_found.remove(faces_to_remove)
				print faces_to_remove, faces_not_found
			except:
				pass
			if faces_not_found == []:
				testEl.interior = True
				break

		return testEl.interior
			
	def print_details(self):

		print "num_elements = %d" % (self.num_elements)
		print "num_surface_elements = %d" % (self.num_surface_elements)
		print "num_interior_elements = %d" % (self.num_interior_elements)
		sleep(1)

		for e in self.element:
			index = self.element.index(e)
			outline = "Element " + str(index) + " "
			if(index < self.num_surface_elements):
				outline += "(Surface): "
			else:
				outline += "(Interior): "
			for n in e.n:
				outline += str(n) + " "

			print outline

	def write_to_file(self, fname):


		print "Writing to " + fname + "..."

		# Write differently depending on format
		base, ext = os.path.splitext(fname)

		if ext == ".vol":
			fout = open(fname, "a")
			fout.write("#  matnr      np      p1      p2      p3      p4\nvolumeelements\n%d\n" % (self.num_elements))
			for e in self.element:
				fout.write("1 %d" % (len(e.n)))
				for n in e.n:
					fout.write(" %d" % (n))
				fout.write("\n")

			fout.write("\n\n")

		else:
			pass

		fout.close()
		print "done!"

	def calc_mass(self, mat, node, scale = 1.0):
	
		mass = 0.0
		index = 0
		for e in self.element:
			mass += e.get_volume(node, scale) * mat.element[index][0]
			index += 1
		return mass

	def reset(self):

		self.element = []
		self.num_elements = 0
		self.num_surface_elements = 0
		self.num_interior_elements = 0
		self.valid = False

class FFEA_element:

	def __init__(self):

		self.reset()

	def set_indices(self, alist):
		
		# Test for correct number of nodes
		if len(alist) != len(self.n):
			print "Incorrect number of nodes for assignment to this element type."
			return

		for i in range(len(alist)):
			self.n[i] = int(alist[i])

	def calc_centroid(self, node):
	
		centroid = np.array([0.0,0.0,0.0])
		for i in self.n:
			centroid += node.pos[i]
			
		return centroid * (1.0 / len(self.n))
	
	def get_linear_face(self, index):
		
		face = FFEA_surface.FFEA_face_tri_lin()
		if index == 0:
			n = [self.n[0], self.n[1], self.n[2]]
		elif index == 1:
			n = [self.n[0], self.n[3], self.n[1]]
		elif index == 2:
			n = [self.n[0], self.n[2], self.n[3]]
		elif index == 3:
			n = [self.n[1], self.n[3], self.n[2]]

		face.set_indices(n)
		return face

	def get_volume(self, node, scale = 1.0):
		e = []
		for i in range(3):
			e.append(node.pos[self.n[i + 1]] - node.pos[self.n[0]])		

		return np.fabs(np.dot(e[2], np.cross(e[1], e[0])) / 6.0) * np.power(scale, 3.0)

	def reset(self):
		
		self.n = []
		self.interior = None

class FFEA_element_tet_lin(FFEA_element):

	def reset(self):

		self.n = [0,1,2,3]
		self.interior = None

class FFEA_element_tet_sec(FFEA_element):

	def reset(self):

		self.n = [0,1,2,3,4,5,6,7,8,9]
		self.interior = None
