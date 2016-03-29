from os import path
from time import sleep

class FFEA_pin:

	def __init__(self, fname):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname):

		print("Loading FFEA pin file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".pin":
			try:
				self.load_pin(fname)
			except:
				print("\tUnable to load FFEA_pin from " + fname + ". Returning empty object...")

		elif ext == ".bsites":
			try:
				self.load_bsites(fname)
			except:
				print("\tUnable to load FFEA_pin from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_pin(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea pinned nodes file" and line != "walrus pinned nodes file":
			print("\tExpected 'ffea pinned nodes file' but found " + line)
			raise TypeError

		num_pinned_nodes = int(fin.readline().split()[1])

		fin.readline()

		# Read pinned nodes now
		pintype = 0	
		while(True):
			line = fin.readline().strip()
			if line == "":
				break
			else:
				self.add_pinned_node(line)

		fin.close()

	def add_pinned_node(self, anint):

		self.index.append(int(anint))
		self.num_pinned_nodes += 1
		
	def print_details(self):

		print "num_pinned_nodes = %d" % (self.num_pinned_nodes)
		sleep(1)

		outline = ""
		for i in self.index:
			outline += str(i) + " "
			
		print outline
	
	def reset(self):

		self.index = []
		self.num_pinned_nodes = 0
