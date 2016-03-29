from os import path
from time import sleep

class FFEA_material:

	def __init__(self, fname = ""):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname):

		print("Loading FFEA material file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".mat":
			try:
				self.load_mat(fname)
			except:
				print("\tUnable to load FFEA_material from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_mat(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea material params file" and line != "walrus material params file":
			print("\tExpected 'ffea material params file' but found " + line)
			raise TypeError

		num_elements = int(fin.readline().split()[1])

		# Read elements now	
		while(True):
			sline = fin.readline().split()
			if len(sline) != 6:
				break
			# Get an element (we want a matrix of values for slicing)
			el = []

			for s in sline:
				el.append(float(s))
			self.add_element(el)

		fin.close()

	def add_element(self, el):

		self.element.append(el)
		self.num_elements += 1

	def get_num_elements(self):

		return len(self.element)

	def print_details(self):

		print "num_elements = %d" % (self.num_elements)

		print "\t\tDensity,Shear Viscosity,Bulk Viscosity,Shear Modulus,Bulk Modulus,Dielectric Constant\n"
		sleep(1)

		index = -1
		for e in self.element:
			index += 1
			outline = "Element " + str(index) + "\t"
			for param in e:
				outline += str(param) + ", "

			print outline

	def reset(self):

		self.element = []
		self.num_elements = 0
