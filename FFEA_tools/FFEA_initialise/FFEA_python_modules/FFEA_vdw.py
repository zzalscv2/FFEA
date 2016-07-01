from os import path
from time import sleep

class FFEA_vdw:

	def __init__(self, fname):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname):

		print("Loading FFEA vdw file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".vdw":
			try:
				self.load_vdw(fname)
				self.valid = True
			except:
				print("\tUnable to load FFEA_vdw from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_vdw(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found. Returning empty object...")
			self.reset()

		# Test format
		line = fin.readline().strip()
		if line != "ffea vdw file" and line != "walrus vdw file":
			print("\tExpected 'ffea vdw file' but found " + line)
			raise TypeError

		num_faces = int(fin.readline().split()[1])

		fin.readline()

		# Read vdw radii now
		while(True):
			line = fin.readline().strip()
			if line == "":
				break
			else:
				self.add_face(line)

		fin.close()

	def set_num_faces(self, num_faces):
		self.num_faces = num_faces
		self.index = [-1 for i in range(num_faces)]
		
	def add_face(self, anint):

		self.index.append(int(anint))
		self.num_faces += 1
		
	def print_details(self):

		print "num_faces = %d" % (self.num_faces)
		sleep(1)

		outline = ""
		for i in self.index:
			outline += "%d " % (i)
			
		print outline
	
	def reset(self):

		self.index = []
		self.num_faces = 0
		self.valid = False
