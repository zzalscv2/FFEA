from os import path

class FFEA_lj:

	def __init__(self, fname = ""):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname):

		print("Loading FFEA lj file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
			return
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".lj":
			try:
				self.load_lj(fname)
				self.valid = True
			except:
				print("\tUnable to load FFEA_lj from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	
	def load_lj(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found. Returning empty object...")
			self.reset()

		# Test format
		line = fin.readline().strip()
		if line != "ffea vdw forcefield params file" and line != "walrus vdw forcefield params file":
			print("\tExpected 'ffea vdw file' but found " + line)
			raise TypeError

		self.num_face_types = int(fin.readline().split()[1])


		# Read interaction matrix now
		for i in range(self.num_face_types):
			sline = fin.readline().split()
			self.interaction.append([])
			for s in sline:
				intline = s.strip()[1:-1].split(",")
				self.interaction[-1].append(FFEA_lj_pair())
				self.interaction[-1][-1].eps = float(intline[0])
				self.interaction[-1][-1].r = float(intline[0])

		fin.close()

	def default(self):

		# Default to steric (ish) params only
		self.interaction = [[FFEA_lj_pair() for i in range(self.num_face_types)] for j in range(self.num_face_types)]
		for i in range(self.num_face_types):
			for j in range(self.num_face_types):
				self.interaction[i][j].eps = 1e-15
				self.interaction[i][j].r = 1e-10

	def write_to_file(self, fname):

		fout = open(fname, "w")
		fout.write("ffea vdw forcefield params file\nnum_vdw_face_types %d\n" % (self.num_face_types))
		for i in self.interaction:
			for j in i:
				fout.write("(%e, %e) " % (j.eps, j.r))
			fout.write("\n")
		fout.close()

	def reset(self):

		self.num_face_types = 7
		self.interaction = []

class FFEA_lj_pair:

	def __init__(self):

		self.reset()

	def reset(self):

		self.r = 0.0
		self.eps = 0.0
