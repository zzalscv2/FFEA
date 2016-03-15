from os import path

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
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".top":
			try:
				load_top(fname)
			except:
				print("\tUnable to load FFEA_topology from " + fname + ". Returning empty object...")

		elif ext == ".vol":
			try:
				load_vol(fname)
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

		self.num_elements = int(fin.readline().split()[1])
		self.num_surface_elements = int(fin.readline().split()[1])
		self.num_interior_elements = int(fin.readline().split()[1])

		fin.readline()

		# Read elements now		
		while(True):
			line = fin.readline()
			try:
				sline = line.split()
			except:
				break

			# Get an element
			el = FFEA_element()
			el.set_indices(sline)
		fin.close()

	def get_num_elements(self):

		return len(self.elements)

	def reset():

		self.elements = []
		self.num_elements = 0
		self.num_surface_elements = 0
		self.num_interior_elements = 0
