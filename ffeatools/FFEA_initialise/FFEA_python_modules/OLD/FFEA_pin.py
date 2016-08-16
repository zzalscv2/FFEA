import numpy as np
import sys, os
import FFEA_node, FFEA_topology

class FFEA_pin:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		if not os.path.exists(fname):
			print "Pin file " + fname  + " not found. Creating empty object....done!"
			self.reset()
			return

		fin = open(fname, "r")

		# Header
		line = fin.readline().rstrip()
		if line != "ffea pinned nodes file" and line != "walrus pinned nodes file":
			print "Error. Expected to read 'ffea pinned nodes file'. This may not be an ffea pinned nodes file"
			return

		# num_pinned_nodes
		try:
			self.num_pinned_nodes = int(fin.readline().split()[1])
			self.node_index = []

		except(ValueError):
			print "Error. Expected to read:"
			print "num_pinned_nodes %d"
			self.reset()
			fin.close()
			return

		# pin indices
		if fin.readline().strip() != "pinned nodes:":
			print "Error. Expected to read 'pinned nodes:' to begin the pin indices section."
			self.reset()
			fin.close()
			return

		for i in range(self.num_pinned_nodes):
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				self.node_index.append(int(line))

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_pinned_nodes = " + str(self.num_pinned_nodes) + "\nnum_pinned_nodes read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected a pin index of the form '%d' for face " + str(i) + ", but found " + line
				self.reset()
				fin.close()
				return
		
		self.node_index = set(self.node_index)
		fin.close()

	def reset(self):
		self.node_index = []
		self.num_pinned_nodes = 0

	def write_to_file(self, fname):

		fout = open(fname, "w")
		fout.write("ffea pinned nodes file\nnum_pinned_nodes %d\npinned nodes:\n" % (self.num_pinned_nodes))
		for i in self.node_index:
			fout.write(str(i) +"\n")

	def get_num_pinned_nodes(self):

		return self.num_pinned_nodes

	def get_num_linear_pinned_nodes(self, top):

		linear_nodes = top.get_linear_nodes()
		return len(linear_nodes & self.node_index)

	def get_node_indices(self):

		return self.node_index

	def pin_radially(self, node, origin, radius):

		self.reset()
		origin = np.array(origin)
		for i in range(len(node.pos)):
			if np.linalg.norm(node.pos[i] - origin) < radius:
				self.num_pinned_nodes += 1
				self.node_index.append(i)

	def pin_linear_radially(self, node, top, origin, radius):

		self.reset()
		origin = np.array(origin)

		# Get list of linear nodes
		linear_index = []
		for e in top.element:
			for i in range(4):
				linear_index.append(e.n[i])

		linear_index = set(linear_index)

		# Now scan over linear indices only
		for i in linear_index:
			if np.linalg.norm(node.pos[i] - origin) < radius:
				self.num_pinned_nodes += 1
				self.node_index.append(i)

	def add_node(self, node_no):

		self.node_index.append(node_no)
		self.num_pinned_nodes += 1
