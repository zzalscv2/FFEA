import FFEA_node, FFEA_surface, FFEA_material, FFEA_topology
import FFEA_vdw, FFEA_pin, FFEA_stokes, FFEA_binding_sites
import numpy as np

class FFEA_blob:

	def __init__(self):
		
		self.reset()

	def load(self, blob_index, conformation_index, nodefname = "", surffname = "", topfname = "", matfname = "", vdwfname = "", pinfname = "", stokesfname = "", bsitesfname = "", motion_state = "DYNAMIC"):

		self.blob_index = blob_index
		self.conformation_index = conformation_index
		self.motion_state = motion_state
		self.surf = FFEA_surface.FFEA_surface(surffname)
		self.node = FFEA_node.FFEA_node(nodefname)

		#if issubclass(self, FFEA_blob):
		self.node.calc_normals(self.surf)

		self.vdw = FFEA_vdw.FFEA_vdw(vdwfname)
		self.bsites = FFEA_binding_sites.FFEA_binding_sites(bsitesfname)

		if self.motion_state == "DYNAMIC":	
			self.top = FFEA_topology.FFEA_topology(topfname)
			self.mat = FFEA_material.FFEA_material(matfname)
			self.pin = FFEA_pin.FFEA_pin(pinfname)
			self.stokes = FFEA_stokes.FFEA_stokes(stokesfname)

	def set_centroid(self, new_cent):
		if new_cent == None:
			return

		cur_cent = self.node.calc_centroid()
		self.node.pos += new_cent - cur_cent

	def apply_rotation(self, rotation):
		if rotation == None:
			return
		
		# Convert to radians
		rotation = np.array(rotation)
		rotation *= np.pi / 180.0
 
		# Build the appropriate rotation matrix
		r = np.array([[0.0 for i in range(3)] for j in range(3)])
		
		# This is an absolute matrix
		if len(rotation) == 9:
			k = 0
			for i in range(3):
				for j in range(3):
					r[i][j] = rotation[k]
					k += 1

		# This is a set of 3 angles x,y,z
		if len(rotation) == 3:
			r[0][0] = np.cos(rotation[1]) * np.cos(rotation[2])
			r[0][1] = np.sin(rotation[0]) * np.sin(rotation[1]) * np.cos(rotation[2]) - np.cos(rotation[0]) * np.sin(rotation[2])
			r[0][2] = np.cos(rotation[0]) * np.sin(rotation[1]) * np.cos(rotation[2]) + np.sin(rotation[0]) * np.sin(rotation[2])
			r[1][0] = np.cos(rotation[1]) * np.sin(rotation[2])
			r[1][1] = np.sin(rotation[0]) * np.sin(rotation[1]) * np.sin(rotation[2]) + np.cos(rotation[0]) * np.cos(rotation[2])
			r[1][2] = np.cos(rotation[0]) * np.sin(rotation[1]) * np.sin(rotation[2]) - np.sin(rotation[0]) * np.cos(rotation[2])
			r[2][0] = -1 * np.sin(rotation[1])
			r[2][1] = np.sin(rotation[0]) * np.cos(rotation[1])
			r[2][2] = np.cos(rotation[0]) * np.cos(rotation[1])

		# Translate the system to the origin, storing original position
		cent = self.node.calc_centroid()
		self.set_centroid(np.array([0.0,0.0,0.0]))
	
		# Apply rotation to each node
		for i in range(self.node.num_nodes):
			self.node.pos[i] = np.dot(r, self.node.pos[i])

		# Translate back
		self.set_centroid(cent)


	def reset(self):

		# FFEA properties
		self.blob_index = 0
		self.conformation_index = 0

		# Structural objects
		self.motion_state = "DYNAMIC"
		self.node = None
		self.surf = None
		self.top = None
		self.mat = None
		self.vdw = None
		self.stokes = None
		self.pin = None
		self.bsites = None
