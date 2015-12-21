import sys, os

class FFEA_springs:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. File " + fname  + " not found. Returning empty FFEA_springs object..."
			return

		# Header
		if fin.readline().rstrip() != "ffea springs file":
			print "Error. Expected to read 'ffea springs file'. This may not be an ffea springs file"
			return

		# num_springs
		try:
			self.num_springs = int(fin.readline().split()[1])

		except(ValueError):
			print "Error. Expected to read:"
			print "num_springs = %d\n"
			self.reset()
			fin.close()
			return

		# springs:
		if fin.readline().strip() != "springs:":
			print "Error. Expected to read 'springs:' to mark the beginning of the springs section\n"
			self.reset()
			fin.close()
			return

		for i in range(self.num_springs):

			try:
				sline = fin.readline().split()
				spring = FFEA_spring()
				spring.k = float(sline[0])
				spring.l = float(sline[1])
				spring.blob_indices = [int(sline[2]), int(sline[3])]
				spring.conformation_indices = [int(sline[4]), int(sline[5])]
				spring.node_indices = [int(sline[6]), int(sline[7])]
			except:
				print "Error. Expected line format:\n"
				print "\tk l blob_index[0], blob_index[1] conformation_index[0] conformation_index[1] node_index[0] node_index[1]\n"
				self.reset()
				fin.close()
				return

			self.spring_array.append(spring)

	def write_to_file(self, fname):

		fout = open(fname, "w")
		fout.write("ffea springs file\nnum_springs %d\nsprings:\n" % (self.num_springs))
		for spring in self.spring_array:
			fout.write("%f %f %d %d %d %d %d %d\n" % (spring.k, spring.l, spring.blob_indices[0], spring.blob_indices[1], spring.conformation_indices[0], spring.conformation_indices[1], spring.node_indices[0], spring.node_indices[1]))

		fout.close()

	def add_spring(self, k, l, bi, ci, ni):

		spring = FFEA_spring()
		spring.k = k
		spring.l = l
		spring.blob_indices = [int(i) for i in bi]
		spring.conformation_indices = [int(i) for i in ci]
		spring.node_indices = [int(i) for i in ni]
	
		self.spring_array.append(spring)
		self.num_springs += 1

	def reset(self):

		self.num_springs = 0
		self.spring_array = []

class FFEA_spring:

	def __init__(self):

		# Initialise stuff
		self.reset()

	def reset(self):

		self.k = 0.0
		self.l = 0.0
		self.blob_indicies = [0,0]
		self.conformation_indices = [0,0]
		self.node_indices = [0,0]
