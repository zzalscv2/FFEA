import numpy as np
import FFEA_node

class FFEA_frame(FFEA_node.FFEA_node):

	# Load with file object already open
	def load_from_traj(self, fo):
		
		while(True):
			prev = fo.tell()
			line = fo.readline().split()
			try:
				self.pos.append([float(line[i]) for i in range(3)])
			except(IndexError):
				break
			except(ValueError):
				fo.seek(prev)
				break

		# Numpy it up for speed
		self.pos = np.array(self.pos)
		self.num_nodes = len(self.pos)
		self.num_surface_nodes = self.num_nodes

	def write_to_traj(self, fo):

		for p in self.pos:
			fo.write("%10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e\n" % (p[0], p[1], p[2], 0, 0, 0, 0, 0, 0, 0))

	# Function to calculate normals at each node, average of connecting faces
	def calc_normals(self, surf):

		self.normal = np.array([[0.0, 0.0, 0.0] for i in range(len(self.pos))])
		for f in surf.face:
			norm = f.get_normal(self)
			for n in f.n:
				self.normal[n] += norm

		for n in self.normal:
			print n
			n /= np.linalg.norm(n)

	def set_step(self, step):
		self.step = step

	def reset(self):
		self.num_nodes = 0
		self.num_surface_nodes = 0
		self.num_interior_nodes = 0
		self.step = 0
		self.pos = []
		self.normal = []
