import sys
import numpy as np

class FFEA_topology:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. Topology file " + fname  + " not found."
			return

		# Header
		if fin.readline().rstrip() != "ffea topology file":
			print "Error. Expected to read 'ffea topology file'. This may not be an ffea topology file"
			return

		# num_elements
		try:
			self.num_elements = int(fin.readline().split()[1])
			self.num_surface_elements = int(fin.readline().split()[1])
			self.num_interior_elements = int(fin.readline().split()[1])

		except(ValueError):
			print "Error. Expected to read:"
			print "num_elements = %d\nnum_surface_elements = %d\nnum_interior_elements = %d"
			self.reset()
			fin.close()
			return			

		# Begin to read elements
		if fin.readline().strip() != "surface elements:":
			print "Error. Expected to read 'surface elements:' to begin the surface elements section."
			self.reset()
			fin.close()
			return

		for i in range(self.num_surface_elements):
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				sline = line.split()

				# Final 6 indices are the secondary element nodes. May need later
				self.element.append(FFEA_element(int(sline[0]), int(sline[1]), int(sline[2]), int(sline[3])))

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_elements = " + str(self.num_elements) + "\nnum_elements read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected a top position of the form '%f %f %f %f' for element " + str(i) + ", but found " + line
				self.reset()
				fin.close()
				return

		if fin.readline().strip() != "interior elements:":
			print "Error. Expected to read 'interior elements:' to begin the interior elements section."
			self.reset()
			fin.close()
			return

		for i in range(self.num_surface_elements, self.num_elements):
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				sline = line.split()

				# Final 6 indices are the secondary element nodes. May need later
				self.element.append(FFEA_element(int(sline[0]), int(sline[1]), int(sline[2]), int(sline[3])))

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_elements = " + str(self.num_elements) + "\nnum_elements read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected a top position of the form '%f %f %f %f' for element " + str(i) + ", but found " + line
				self.reset()
				fin.close()
				return

	def reset(self):
		self.num_elements = 0
		self.element = []

	def get_linear_nodes(self):
		
		linear_nodes = []
		for elem in self.element:
			for i in elem.n:
				linear_nodes.append(i)

		return set(linear_nodes)

	def get_num_linear_nodes(self):
		
		return len(self.get_linear_nodes())

class FFEA_element:

	def __init__(self, n0, n1, n2, n3):
	
		self.n = [n0, n1, n2, n3]
		self.normal = []
		self.face_centroid = []

	def calc_centroid(self, ffea_node):
		
		centroid = np.array([0.0,0.0,0.0])
		for index in self.n:
			centroid += ffea_node.pos[index]

		return centroid * (1.0/4.0)

	def calc_normals(self, node):

		self.normal = [self.calc_face_normal(i, node) for i in range(4)]

	def calc_face_centroids(self, node):

		self.face_centroid = [self.calc_face_centroid(i, node) for i in range(4)]

	def calc_face_normal(self, index, node):

		if index == 0:

			# n[0], n[1], n[2]
			v0 = node.pos[self.n[1]] - node.pos[self.n[0]]
			v1 = node.pos[self.n[2]] - node.pos[self.n[0]]

		elif index == 1:
			
			# n[0], n[1], n[3]
			v0 = node.pos[self.n[3]] - node.pos[self.n[0]]
			v1 = node.pos[self.n[1]] - node.pos[self.n[0]]

		elif index == 2:
			
			# n[0], n[2], n[3]
			v0 = node.pos[self.n[2]] - node.pos[self.n[0]]
			v1 = node.pos[self.n[3]] - node.pos[self.n[0]]

		elif index == 3:
			
			# n[1], n[2], n[3]
			v0 = node.pos[self.n[1]] - node.pos[self.n[2]]
			v1 = node.pos[self.n[3]] - node.pos[self.n[2]]

		v = np.cross(v0, v1)
		return v / np.linalg.norm(v)

	def calc_face_centroid(self, index, node):

		if index == 0:

			# n[0], n[1], n[2]
			return (node.pos[self.n[0]] + node.pos[self.n[1]] + node.pos[self.n[2]]) / 3.0

		elif index == 1:

			# n[0], n[1], n[3]
			return (node.pos[self.n[0]] + node.pos[self.n[1]] + node.pos[self.n[3]]) / 3.0

		elif index == 2:

			# n[0], n[2], n[3]
			return (node.pos[self.n[0]] + node.pos[self.n[2]] + node.pos[self.n[3]]) / 3.0

		elif index == 3:

			# n[1], n[2], n[3]
			return (node.pos[self.n[1]] + node.pos[self.n[2]] + node.pos[self.n[3]]) / 3.0

	def contains(self, pos, node):
		
		# And separations of point from centroid
		sep = [i - pos for i in self.face_centroid]
		sep_norm = [s / np.linalg.norm(s) for s in sep]

		# Now, the check
		for i in range(4):
			if np.dot(sep_norm[i], self.normal[i]) < 0:
				return False

		return True
