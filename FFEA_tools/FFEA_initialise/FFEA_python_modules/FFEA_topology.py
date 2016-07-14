import os
from time import sleep
import numpy as np

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
			self.add_element(el, eltype = eltype)

		fin.close()

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
		num_interior_elements = num_elements
		num_surface_elements = 0

		# Read elements now
		eltype = 1
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
			self.add_element(el, eltype = eltype)

		fin.close()

	def add_element(self, el, eltype = 0):

		self.element.append(el)
		self.num_elements += 1

		if eltype == 0:
			self.num_surface_elements += 1
		else:
			self.num_interior_elements += 1

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
		
	def reset(self):
		
		self.n = []

class FFEA_element_tet_lin(FFEA_element):

	def reset(self):

		self.n = [0,1,2,3]

class FFEA_element_tet_sec(FFEA_element):

	def reset(self):

		self.n = [0,1,2,3,4,5,6,7,8,9]
