from os import path
from time import sleep

class FFEA_stokes:

	def __init__(self, fname):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname):

		print("Loading FFEA stokes file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".stokes":
			try:
				self.load_stokes(fname)
			except:
				print("\tUnable to load FFEA_stokes from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_stokes(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea stokes radii file" and line != "walrus stokes radii file":
			print("\tExpected 'ffea stokes radii file' but found " + line)
			raise TypeError

		num_nodes = int(fin.readline().split()[1])

		fin.readline()

		# Read stokes radii now
		while(True):
			line = fin.readline().strip()
			if line == "":
				break
			else:
				self.add_node(line)

		fin.close()

	def add_node(self, afloat):

		self.radius.append(float(afloat))
		self.num_nodes += 1
		
	def print_details(self):

		print "num_nodes = %d" % (self.num_nodes)
		sleep(1)

		outline = ""
		for rad in self.radius:
			outline += "%6.3f\n" % (rad)
			
		print outline
	
	def reset(self):

		self.radius = []
		self.num_nodes = 0
