from os import path
from time import sleep

class FFEA_springs:

	def __init__(self, fname = ""):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname):

		print("Loading FFEA springs file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".springs":
			try:
				self.load_springs(fname)
			except:
				print("\tUnable to load FFEA_springs from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_springs(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea springs file" and line != "walrus springs file":
			print("\tExpected 'ffea springs file' but found " + line)
			raise TypeError

		num_springs = int(fin.readline().split()[1])

		fin.readline()

		# Read springs now	
		while(True):
			sline = fin.readline().split()

			# Get a spring
			spring = FFEA_spring()
			try:
				spring.set_properties(sline[0], sline[1], sline[2:4], sline[4:6], sline[6:8])

			except(IndexError):
				break

			except:
				break

			self.add_spring(spring)

		fin.close()

	def add_spring(self, spring):

		self.spring.append(spring)
		self.num_springs += 1

	def get_num_springs(self):

		return self.num_springs

	def print_details(self):

		print "num_springs = %d" % (self.num_springs)
		sleep(1)

		print "\n\t k\t\tl\t\tblob\tconf\tnode"
		for s in self.spring:
			index = self.spring.index(s)
			outline = "Spring " + str(index) + " "
			outline += "%e\t%e\t%d %d\t%d %d\t%d %d" % (s.k, s.l, s.blob_index[0], s.blob_index[1], s.conformation_index[0], s.conformation_index[1], s.node_index[0], s.node_index[1])

			print outline

	def reset(self):

		self.spring = []
		self.num_springs = 0

class FFEA_spring:

	def __init__(self):

		self.reset()

	def set_properties(self, k, l, bin, cin, nin):

		try:
			self.k = float(k)
			self.l = float(l)
			self.blob_index = [int(i) for i in bin]
			self.conformation_index = [int(i) for i in cin]
			self.node_index = [int(i) for i in nin]
		except:
			raise

	def reset(self):
		
		self.k = 0
		self.l = 0
		self.blob_index = []
		self.conformation_index = []
		self.node_index = []
