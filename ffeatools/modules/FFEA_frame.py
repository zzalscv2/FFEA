import numpy as np
import FFEA_node

class FFEA_frame(FFEA_node.FFEA_node):

	# Load with file object already open
	def load_from_traj(self, fo):
		
		start = fo.tell()
		
		while(True):
			prev = fo.tell()
			line = fo.readline().split()
			try:
				self.pos.append([float(line[i]) for i in range(3)])
				
			except(IndexError):
			
				# EOF
				fo.seek(start)
				return 1

			except(ValueError):
				if line[0] == "*" or line[0] == "Blob":
				
					# Ready for next frame
					fo.seek(prev)
					break
				else:
					# Halfway through a written frame
					fo.seek(start)
					return 1
					
		# Numpy it up for speed
		self.pos = np.array(self.pos)
		self.num_nodes = len(self.pos)
		self.num_surface_nodes = self.num_nodes
		
		return 0
		
	def build_from_node(self, node):
	
		self.num_nodes = node.num_nodes
		self.num_surface_nodes = node.num_surface_nodes
		self.num_interior_nodes = node.num_interior_nodes
		self.pos = node.pos
		
	def write_to_traj(self, fo):

		for p in self.pos:
			fo.write("%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n" % (p[0], p[1], p[2], 0, 0, 0, 0, 0, 0, 0))

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
