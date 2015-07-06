import sys, os

class FFEA_pin:

	def __init__(self, pin_fname):
		
		# Open file and check initial stuff
		try:
			fin = open(pin_fname, "r")
		except(IOError):
			self.num_pinned_nodes = 0
			self.pinned_node_index = []
			return

		line = fin.readline().strip()
		if line != "ffea pinned nodes file" and line != "walrus pinned nodes file":
			sys.exit("Error. Expected but didn't find 'ffea pinned nodes file'\n")
		
		self.num_pinned_nodes = int(fin.readline().split()[1].strip())
		self.pinned_node_index = [0 for i in range(self.num_pinned_nodes)]
		fin.readline()	#'pinned nodes:'
		
		# Read pinned node indices
		for i in range(self.num_pinned_nodes):
			self.pinned_node_index[i] = int(fin.readline().strip())

		fin.close()

	def reset(self):
		self.num_pinned_nodes = 0
		self.pinned_node_index = []

	def add_node(self, index):
		self.num_pinned_nodes += 1
		self.pinned_node_index.append(index)

	def write_to_file(self, fname):

		fout = open(fname, "w")
		fout.write("ffea pinned nodes file\n")
		fout.write("num_pinned_nodes %d\n" % (self.num_pinned_nodes))
		fout.write("pinned nodes:\n")
		for i in range(self.num_pinned_nodes):
			fout.write("%d\n" % (self.pinned_node_index[i]))

		fout.close()