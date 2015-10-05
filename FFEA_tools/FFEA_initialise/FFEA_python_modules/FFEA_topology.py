import sys

class FFEA_topology:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. File " + fname  + " not found."
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
