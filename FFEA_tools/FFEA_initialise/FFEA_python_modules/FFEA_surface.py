import sys
import numpy as np

class FFEA_surface:

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
		if fin.readline().rstrip() != "ffea surface file":
			print "Error. Expected to read 'ffea surface file'. This may not be an ffea surf file"
			return

		# num_faces
		try:
			self.num_faces = int(fin.readline().split()[1])

		except(ValueError):
			print "Error. Expected to read:"
			print "num_surface_faces = %d"
			self.reset()
			fin.close()
			return			

		# Begin to read faces
		if fin.readline().strip() != "faces:":
			print "Error. Expected to read 'faces:' to begin the surface section."
			self.reset()
			fin.close()
			return

		for i in range(self.num_faces):
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				sline = line.split()

				# First index is the containing element. May need later
				self.face.append(FFEA_face(int(sline[1]), int(sline[2]), int(sline[3])))

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_faces = " + str(self.num_faces) + "\nnum_faces read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected a surf position of the form '%f %f %f %f' for face " + str(i) + ", but found " + line
				self.reset()
				fin.close()
				return

	def calculate_structure_dimensions(self, node):

		xlim = [float("inf"), -1 * float("inf")]
		ylim = [float("inf"), -1 * float("inf")]
		zlim = [float("inf"), -1 * float("inf")]

		for f in self.face:
			for n in f.n:
				if node.pos[n][0] < xlim[0]:
					xlim[0] = node.pos[n][0]
				elif node.pos[n][0] > xlim[1]:
					xlim[1] = node.pos[n][0]

				if node.pos[n][1] < ylim[0]:
					ylim[0] = node.pos[n][1]
				elif node.pos[n][1] > ylim[1]:
					ylim[1] = node.pos[n][1]

				if node.pos[n][2] < zlim[0]:
					zlim[0] = node.pos[n][2]
				elif node.pos[n][2] > zlim[1]:
					zlim[1] = node.pos[n][2]

		return np.array([xlim[1]-xlim[0], ylim[1]-ylim[0], zlim[1]-zlim[0]]), np.array([xlim, ylim, zlim])

	def reset(self):
		self.num_faces = 0
		self.face = []

class FFEA_face:

	def __init__(self, n0, n1, n2):
	
		self.n = [n0, n1, n2]

	def calc_centroid(self, node):

		centroid = np.array([0.0,0.0,0.0])
		for n in self.n:
			centroid += node.pos[n]
		
		return centroid * 1.0/3.0

	def calc_area(self, node):
		
		v1 = node.pos[self.n[1]] - node.pos[self.n[0]]
		v2 = node.pos[self.n[2]] - node.pos[self.n[0]]
		return 0.5 * np.linalg.norm(np.cross(v1,v2))

	def get_normal(self, node):

		v1 = node.pos[self.n[2]] - node.pos[self.n[1]]
		v2 = node.pos[self.n[1]] - node.pos[self.n[0]]
		norm = np.cross(v1,v2)
		return norm * 1.0 / np.linalg.norm(norm)
