from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import math, os
import numpy as np
import FFEA_material

class Blob:
	def __init__(self, energy_thresh=1.0e6):
		self.num_elements = 0
		self.offset = [0.0, 0.0, 0.0]
		self.init_centroid = None
		self.init_rot = None
		self.topology = []
		self.no_topology = False
		self.num_nodes = 0
		self.num_surface_faces = 0
		self.num_binding_sites = 0
		self.active_binding_site = -1
		self.num_frames = 0
		self.frames = []
		self.num_frames = 0
		self.state = "DYNAMIC"
		self.id_num = -1
		self.hide_blob = False
		self.num_pinned_nodes = 0
		self.pinned_nodes = []
		self.min_length = None
		self.linear_nodes_only = []
		self.calculated_linear_nodes = False
		self.calculated_first_frame_volumes = False
		self.calculated_first_frame_J_inv = False
		self.do_Fij = False
		self.energy_thresh = energy_thresh
		self.scale = 1.0
		self.global_scale = 1.0
		self.normalcolor = [32 / 255.0, 178 / 255.0, 170 / 255.0]

	def load(self, idnum, blob_index, conformation_index, nodes_fname, top_fname, surf_fname, vdw_fname, scale, blob_state, blob_pinned, blob_mat, binding_fname, blob_centroid_pos, blob_rotation, ffea_fpath = "."):
		self.id_num = idnum
		self.blob_index = blob_index
		self.conformation_index = conformation_index
                if os.path.isabs(nodes_fname) == False:
                     self.nodes_fname = os.path.join(ffea_fpath, nodes_fname)
                else:
		     self.nodes_fname = nodes_fname
		self.scale = scale
                self.ffea_fpath = ffea_fpath
		
		if blob_centroid_pos != None:
			self.init_centroid = blob_centroid_pos
			self.offset = blob_centroid_pos
			
		if blob_rotation != None:
			self.init_rot = blob_rotation

		self.state = blob_state

		if top_fname == "":
			print "No topology file provided."
			self.no_topology = True
		else:
                        if os.path.isabs(top_fname) == False:
                             top_fname = os.path.join(ffea_fpath, top_fname)

			self.load_topology(top_fname)
			self.no_topology = False

                if os.path.isabs(surf_fname) == False:
                     surf_fname = os.path.join(ffea_fpath, surf_fname)
		self.load_surface(surf_fname)

		if self.state == "DYNAMIC":
			self.mat = FFEA_material.FFEA_material(blob_mat)
		else:
			self.mat = None

		if vdw_fname == "":
			print "No vdw file provided. Creating a zero array."
			self.vdw = [-1 for i in xrange(self.num_surface_faces)]
		else:
                        if os.path.isabs(vdw_fname) == False:
                             vdw_fname = os.path.join(ffea_fpath, vdw_fname)
			try:
				self.load_vdw(vdw_fname)
			except(IOError):
				print "Vdw file '" + vdw_fname + "' not found. Creating a zero array."
				self.vdw = [-1 for i in xrange(self.num_surface_faces)]
			except:
				print "Unknow error in 'self.load_vdw()'. Creating a zero array."
				self.vdw = [-1 for i in xrange(self.num_surface_faces)]

		self.hidden_face = [-1 for i in xrange(self.num_surface_faces)]
		for i in xrange(self.num_surface_faces):
			if self.vdw[i] == -2:
				self.hidden_face[i] = 1

		if blob_pinned == "":
			print "No pinned nodes file provided."
			self.num_pinned_nodes = 0
		else:
                        if os.path.isabs(blob_pinned) == False:
                             blob_pinned = os.path.join(ffea_fpath, blob_pinned)

			self.load_pinned_nodes(blob_pinned)

		if binding_fname == "":
			print "No binding sites will be loaded."
		
		else:
                        if os.path.isabs(binding_fname) == False:
                             binding_fname = os.path.join(ffea_fpath, binding_fname)
			self.load_binding_sites(binding_fname)

	def load_topology(self, top_fname):
		print "Reading in topology file " + top_fname
       
		top = open(top_fname, "r")
		line = top.readline().rstrip()
		if line != "ffea topology file" and line != "walrus topology file":
			print "Error: Topology file " + top_fname + " missing 'ffea topology file' first line"
			return

		line = top.readline().split()
		self.num_elements = int(line[1])
		print "num_elements = ", self.num_elements

		line = top.readline().split()
		num_surface_elements = int(line[1])
		print "num_surface_elements = ", num_surface_elements

		line = top.readline().split()
		num_interior_elements = int(line[1])
		print "num_interior_elements = ", num_interior_elements

		line = top.readline().rstrip()
		if line != "surface elements:":
			print "Error: Topology file " + top_fname + " missing 'surface elements:' line"
			return

		for n in xrange(num_surface_elements):
			line = top.readline().split()
			self.topology.append([int(line[0]), int(line[1]), int(line[2]), int(line[3]), int(line[4]), int(line[5]), int(line[6]), int(line[7]), int(line[8]), int(line[9])])

		line = top.readline().rstrip()
		if line != "interior elements:":
			print "Error: Topology file " + top_fname + " missing 'interior elements:' line"
			return

		for n in xrange(num_interior_elements):
			line = top.readline().split()
			self.topology.append([int(line[0]), int(line[1]), int(line[2]), int(line[3]), int(line[4]), int(line[5]), int(line[6]), int(line[7]), int(line[8]), int(line[9])])

		top.close()
		print "Finished reading in topology file " + top_fname

	def load_surface(self, surf_fname):
		print "Reading in surface file " + surf_fname
		surf = open(surf_fname, "r")
		line = surf.readline().rstrip()
		if line != "ffea surface file" and line != "walrus surface file":
			print "Error: surface file " + surf_fname + " missing 'ffea surface file' first line"
			return

		line = surf.readline().split()
		self.num_surface_faces = int(line[1])
		print "num_surface_faces = ", self.num_surface_faces

		line = surf.readline().rstrip()
		if line != "faces:":
			print "Error: surface file " + surf_fname + " missing 'faces:' line"
			return

		self.surface = []
		for n in xrange(self.num_surface_faces):
			line = surf.readline().split()
			self.surface.append([int(line[0]), int(line[1]), int(line[2]), int(line[3])])

		surf.close()
		print "Finished reading in surface file " + surf_fname

	def load_pinned_nodes(self, pin_fname):
		print "Reading in pinned nodes file " + pin_fname
		pin = open(pin_fname, "r")
		line = pin.readline().rstrip()
		if line != "ffea pinned nodes file" and line != "walrus pinned nodes file":
			print "Error: pinned nodes file " + pin_fname + " missing 'ffea pinned nodes file' first line"
			return

		line = pin.readline().split()
		self.num_pinned_nodes = int(line[1])
		print "num_pinned_nodes = ", self.num_pinned_nodes

		line = pin.readline().rstrip()
		if line != "pinned nodes:":
			print "Error: pinned nodes file " + pin_fname + " missing 'pinned nodes:' line"
			return

		self.pinned_nodes = []
		for n in xrange(self.num_pinned_nodes):
			line = pin.readline()
			self.pinned_nodes.append(int(line))

		pin.close()
		print "Finished reading in pinned nodes file " + pin_fname

	def load_vdw(self, vdw_fname):
		print "Reading in vdw file " + vdw_fname
		vdw_file = open(vdw_fname, "r")
		line = vdw_file.readline().rstrip()
		if line != "ffea vdw file" and line != "walrus vdw file":
			print "Error: vdw file " + vdw_fname + " missing 'ffea vdw file' first line"
			return

		line = vdw_file.readline().split()
		num_vdw_faces = int(line[1])
		print "num_faces according to vdw file = ", num_vdw_faces

		if num_vdw_faces != self.num_surface_faces:
			print "Error. Number of faces in vdw file (" + str(num_vdw_faces) + ") does not match number of faces in surface file (" + str(self.num_surface_faces) + ")"
			return

		line = vdw_file.readline().rstrip()
		if line != "vdw params:":
			print "Error: vdw file " + vdw_fname + " missing 'vdw params:' line"
			return

		self.vdw = []
		for n in xrange(self.num_surface_faces):
			line = vdw_file.readline()
			self.vdw.append(int(line))

		vdw_file.close()
		print "Finished reading in vdw file " + vdw_fname

	def load_binding_sites(self, fname):

		print "Reading in binding sites file " + fname
		fin = open(fname, "r")
		line = fin.readline().strip()
		if line != "ffea binding sites file":
			print "Error: binding site file " + fname + " missing 'ffea binding site file' first line"
			return

		self.num_binding_sites = int(fin.readline().split()[1])
		print "num_binding_sites according to binding site file = ", self.num_binding_sites

		fin.readline()		
		self.binding_site = []
		self.binding_site_type = [[-1, -1] for i in range(self.num_surface_faces)]
		
		for i in range(self.num_binding_sites):
			asite = []

			# Type and size first
			sline = fin.readline().split()
			site_type = int(sline[1])
			num_faces = int(sline[3])

			# Now get faces
			sline = fin.readline().strip().split()[1:]
			if(len(sline) != num_faces):
				sys.exit("Error. Specified 'num_faces' not equal to num_faces provided.")
					
			for j in range(len(sline)):
				
				# Ignore num_faces (can't be arsed with a new class)
				asite.append(int(sline[j]))
				
				# Stores type and index
				self.binding_site_type[asite[-1]][1] = i

			self.binding_site.append(asite)

		fin.close()
		print "Finished reading in binding site file " + fname

	def write_vdw(self, vdw_fname):
		print "Writing vdw file " + vdw_fname
		vdw_file = open(vdw_fname, "w")
		vdw_file.write("ffea vdw file\n")
		vdw_file.write("num_faces " + str(self.num_surface_faces) + "\n")
		vdw_file.write("vdw params:\n")
		for n in xrange(self.num_surface_faces):
			vdw_file.write(str(self.vdw[n]) + "\n")
		vdw_file.close()
		print "Finished writing vdw file " + vdw_fname

	def load_frame(self, traj_file):

		# Inactive conf
		if traj_file == None:
			self.frames.append(None)
			self.num_frames += 1
			return

		# There is a trajectory! Ignore offset
		self.offset = [0.0, 0.0, 0.0]

		# skip blob, step line
		line = traj_file.readline()
		blob_state = traj_file.readline().rstrip()
		
		nodes = []
		centroid_x = 0.0
		centroid_y = 0.0
		centroid_z = 0.0
		#cdef int n
		for n in xrange(self.num_nodes):
			line = traj_file.readline().split()
			el_nodes = [float(line[i])* self.global_scale for i in xrange(10)]
			nodes.append(el_nodes)
		
			centroid_x += el_nodes[0]
			centroid_y += el_nodes[1]
			centroid_z += el_nodes[2]

		centroid_x /= self.num_nodes
		centroid_y /= self.num_nodes
		centroid_z /= self.num_nodes
		
		# Calculate average normal at each node (for gl lighting effects)
		normal_list = [[0.0, 0.0, 0.0] for i in xrange(self.num_nodes)]
		#cdef int f
		for f in xrange(self.num_surface_faces):
			# get node indices of this face
			i1 = self.surface[f][1]
			i2 = self.surface[f][2]
			i3 = self.surface[f][3]
		
			# get the normal of the face
			norm = self.calc_normal(nodes[i1], nodes[i2], nodes[i3])
		
			normal_list[i1][0] += norm[0]
			normal_list[i1][1] += norm[1]
			normal_list[i1][2] += norm[2]
		
			normal_list[i2][0] += norm[0]
			normal_list[i2][1] += norm[1]
			normal_list[i2][2] += norm[2]
		
			normal_list[i3][0] += norm[0]
			normal_list[i3][1] += norm[1]
			normal_list[i3][2] += norm[2]
		
		
		self.frames.append(Frame(blob_state, nodes, normal_list, centroid_x, centroid_y, centroid_z))
		self.num_frames += 1

		if self.calculated_linear_nodes == False:
			if self.no_topology == False:
				self.linear_nodes_only = []
				for i in range(self.num_nodes):
					for el in self.topology:
						if i in el[0:4]:
							self.linear_nodes_only.append(i)
							break
				print "Found", len(self.linear_nodes_only), "linear nodes."
				self.calculated_linear_nodes = True	

		if self.calculated_first_frame_volumes == False:
			self.first_frame_vol_list = []
			if self.no_topology == False:
				print "Calculating volumes of elements in first frame..."

				# Find first frame available
				for i in range(len(self.frames)):
					if self.frames[i] != None:
						frame_index = i
						break

				for el in xrange(self.num_elements):
					# Get the indices of the 4 nodes of this tetrahedron
					i1 = self.topology[el][0]
					i2 = self.topology[el][1]
					i3 = self.topology[el][2]
					i4 = self.topology[el][3]

					# Get the nodes
					n1a = self.frames[frame_index].node_list[i1]
					n2a = self.frames[frame_index].node_list[i2]
					n3a = self.frames[frame_index].node_list[i3]
					n4a = self.frames[frame_index].node_list[i4]
					n1 = n1a[0:3]
					n2 = n2a[0:3]
					n3 = n3a[0:3]
					n4 = n4a[0:3]
	
					# get the volume of the elements
					self.first_frame_vol_list.append(self.get_element_volume(n1, n2, n3, n4))
				print "...Done."
				self.calculated_first_frame_volumes = True


		if self.calculated_first_frame_J_inv == False:
			self.first_frame_J_inv = []
			if self.no_topology == False:
				print "Calculating J inv of elements in first frame..."

				# Find first frame available
				for i in range(len(self.frames)):
					if self.frames[i] != None:
						frame_index = i
						break

				for el in xrange(self.num_elements):
					# Get the indices of the 4 nodes of this tetrahedron
					i1 = self.topology[el][0]
					i2 = self.topology[el][1]
					i3 = self.topology[el][2]
					i4 = self.topology[el][3]

					# Get the nodes
					n1a = self.frames[frame_index].node_list[i1]
					n2a = self.frames[frame_index].node_list[i2]
					n3a = self.frames[frame_index].node_list[i3]
					n4a = self.frames[frame_index].node_list[i4]
					n1 = n1a[0:3]
					n2 = n2a[0:3]
					n3 = n3a[0:3]
					n4 = n4a[0:3]
	
					# get the J_inv of the elements
					self.first_frame_J_inv.append(self.get_J_inv(n1, n2, n3, n4))
				print "...Done."
				self.calculated_first_frame_J_inv = True

	def delete_all_frames(self):
		self.num_frames = 0
		self.frames = []

	def load_nodes_file_as_frame(self):
		print "Reading in nodes file " + self.nodes_fname
		nodes_file = open(self.nodes_fname, "r")
		line = nodes_file.readline().rstrip()
		if line != "ffea node file" and line != "walrus node file":
			print "Error: Nodes file " + self.nodes_fname + " missing 'ffea node file' first line"
			return

		line = nodes_file.readline().split()
		self.num_nodes = int(line[1])
		print "num_nodes = ", self.num_nodes

		line = nodes_file.readline().split()
		num_surface_nodes = int(line[1])
		print "num_surface_nodes = ", num_surface_nodes

		line = nodes_file.readline().split()
		num_interior_nodes = int(line[1])
		print "num_interior_nodes = ", num_interior_nodes

		line = nodes_file.readline().rstrip()
		if line != "surface nodes:":
			print "Error: nodes file " + self.nodes_fname + " missing 'surface nodes:' line"
			print line
			return

		print "Scaling by " + str(self.scale * self.global_scale)

		nodes = []
		centroid_x = 0.0
		centroid_y = 0.0
		centroid_z = 0.0
		for n in xrange(self.num_nodes):
			line = nodes_file.readline().split()
			if "interior" in line[0]:
				print "Skipping 'interior nodes:' line"
				line = nodes_file.readline().split()
			
			el_nodes = [float(line[i])*self.scale*self.global_scale for i in xrange(3)]
			el_nodes.extend([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
			nodes.append(el_nodes)
		
			centroid_x += el_nodes[0]
			centroid_y += el_nodes[1]
			centroid_z += el_nodes[2]

		centroid_x/= self.num_nodes
		centroid_y/= self.num_nodes
		centroid_z/= self.num_nodes

		nodes_file.close()
		print "Finished reading in nodes file " + self.nodes_fname

		if self.init_centroid != None:
			print "Moving to starting position..."
			translate = [self.init_centroid[0]*self.scale*self.global_scale - centroid_x, self.init_centroid[1]*self.scale*self.global_scale - centroid_y, self.init_centroid[2]*self.scale*self.global_scale - centroid_z]
					
			for i in range(len(nodes)):
				for j in range(3):
					nodes[i][j] += translate[j]

			centroid_x = self.init_centroid[0]*self.scale*self.global_scale
			centroid_y = self.init_centroid[1]*self.scale*self.global_scale
			centroid_z = self.init_centroid[2]*self.scale*self.global_scale
			print "...done!\n"
		
		# Apply a rotation if necessary
		if self.init_rot != None:
			
			# Move centroid to origin
			for i in range(len(nodes)):
				nodes[i][0] -= centroid_x
				nodes[i][1] -= centroid_y
				nodes[i][2] -= centroid_z
				
			# Rotate and move back to not the origin!
			if len(self.init_rot) == 3:
				for i in range(3):
					self.init_rot[i] = np.radians(self.init_rot[i])
	
				R = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
				R[0][0] = np.cos(self.init_rot[1]) * np.cos(self.init_rot[2])
				R[0][1] = np.sin(self.init_rot[0]) * np.sin(self.init_rot[1]) * np.cos(self.init_rot[2]) - np.cos(self.init_rot[0]) * np.sin(self.init_rot[2])
				R[0][2] = np.cos(self.init_rot[0]) * np.sin(self.init_rot[1]) * np.cos(self.init_rot[2]) + np.sin(self.init_rot[0]) * np.sin(self.init_rot[2])
				R[1][0] = np.cos(self.init_rot[1]) * np.sin(self.init_rot[2])
				R[1][1] = np.sin(self.init_rot[0]) * np.sin(self.init_rot[1]) * np.sin(self.init_rot[2]) + np.cos(self.init_rot[0]) * np.cos(self.init_rot[2])
				R[1][2] = np.cos(self.init_rot[0]) * np.sin(self.init_rot[1]) * np.sin(self.init_rot[2]) - np.sin(self.init_rot[0]) * np.cos(self.init_rot[2])
				R[2][0] = -1 * np.sin(self.init_rot[1])
				R[2][1] = np.sin(self.init_rot[0]) * np.cos(self.init_rot[1])
				R[2][2] = np.cos(self.init_rot[0]) * np.cos(self.init_rot[1])

				for i in range(len(nodes)):
					x = nodes[i][0]
					y = nodes[i][1]
					z = nodes[i][2]

					nodes[i][0] = x * R[0][0] + y * R[0][1] + z * R[0][2] + self.offset[0]
					nodes[i][1] = x * R[1][0] + y * R[1][1] + z * R[1][2] + self.offset[1]
					nodes[i][2] = x * R[2][0] + y * R[2][1] + z * R[2][2] + self.offset[2]

			elif len(self.init_rot) == 9:

				# Explicit matrix
				for i in range(len(nodes)):
					x = nodes[i][0]
					y = nodes[i][1]
					z = nodes[i][2]

					nodes[i][0] = x * self.init_rot[0] + y * self.init_rot[1] + z * self.init_rot[2] + self.offset[0]
					nodes[i][1] = x * self.init_rot[3] + y * self.init_rot[4] + z * self.init_rot[5] + self.offset[1]
					nodes[i][2] = x * self.init_rot[6] + y * self.init_rot[7] + z * self.init_rot[8] + self.offset[2]

			print "...done!\n"
				

		# Calculate average normal at each node (for gl lighting effects)
		print "Calculating node normals for lighting..."
		normal_list = [[0.0, 0.0, 0.0] for i in xrange(self.num_nodes)]
		for f in xrange(self.num_surface_faces):
			# get node indices of this face
			i1 = self.surface[f][1]
			i2 = self.surface[f][2]
			i3 = self.surface[f][3]
		
			# get the normal of the face
			norm = self.calc_normal(nodes[i1], nodes[i2], nodes[i3])
		
			normal_list[i1][0] += norm[0]
			normal_list[i1][1] += norm[1]
			normal_list[i1][2] += norm[2]
		
			normal_list[i2][0] += norm[0]
			normal_list[i2][1] += norm[1]
			normal_list[i2][2] += norm[2]
		
			normal_list[i3][0] += norm[0]
			normal_list[i3][1] += norm[1]
			normal_list[i3][2] += norm[2]
		print "Done."
		
		self.frames.append(Frame(self.get_state(), nodes, normal_list, centroid_x, centroid_y, centroid_z))
		self.num_frames += 1

		if self.no_topology == False:
			print "Calculating surface"
			self.linear_nodes_only = []
			for i in range(self.num_nodes):
				for el in self.topology:
					if i in el[0:4]:
						self.linear_nodes_only.append(i)
						break
			print "Found", len(self.linear_nodes_only), "linear nodes."

		if self.calculated_linear_nodes == False:
			if self.no_topology == False:
				self.linear_nodes_only = []
				for i in range(self.num_nodes):
					for el in self.topology:
						if i in el[0:4]:
							self.linear_nodes_only.append(i)
							break
				print "Found", len(self.linear_nodes_only), "linear nodes."
				self.calculated_linear_nodes = True
#			else:
#				self.linear_nodes_only = []
#				for i in range(self.num_nodes):
#					for f in self.surface:
#						if i in f[0:3]:
#							self.linear_nodes_only.append(i)
#							break
#				print "Found", len(self.linear_nodes_only), "linear nodes."
#				self.calculated_linear_nodes = True

	def set_scale(self, scale):
		self.scale = scale

	def set_global_scale(self, global_scale):
		self.global_scale = global_scale

	def calc_normal(self, n1, n2, n3):
		ax = n2[0] - n1[0]
		ay = n2[1] - n1[1]
		az = n2[2] - n1[2]
		bx = n3[0] - n1[0]
		by = n3[1] - n1[1]
		bz = n3[2] - n1[2]

		return [az * by - ay * bz, ax * bz - az * bx, ay * bx - ax * by]

	def get_centroid(self, i):

		if self.state == "STATIC":
			i = 0

		if self.num_frames == 0:
			return 0.0, 0.0, 0.0

		if self.frames[i] == None:
			return None, None, None

		if i < 0:
			i == 0

		elif i >= self.num_frames:
			i = self.num_frames - 1

		f = self.frames[i]
		return f.centroid_x, f.centroid_y, f.centroid_z

	def draw_frame(self, i, display_flags):

		if self.state == "STATIC":
			i = 0

		if self.num_frames == 0:
			return

		if self.frames[i] == None:
			return

		if self.hide_blob == True:
			return

		if i < 0:
			i == 0
		elif i >= self.num_frames:
			i = self.num_frames - 1

		if display_flags['hide_frozen'] == 1:
			if self.get_frame_state(i) == "FROZEN":
				return
		
		#cdef int f

		if display_flags['vdw_edit_mode'] == 1 and self.id_num == display_flags['selected_index']:
			glBegin(GL_TRIANGLES)
			for f in range(self.num_surface_faces):

				n1 = self.frames[i].node_list[self.surface[f][1]][0:3]
				n2 = self.frames[i].node_list[self.surface[f][2]][0:3]
				n3 = self.frames[i].node_list[self.surface[f][3]][0:3]

				norm1 = self.frames[i].normal_list[self.surface[f][1]]
				norm2 = self.frames[i].normal_list[self.surface[f][2]]
				norm3 = self.frames[i].normal_list[self.surface[f][3]]

				if self.vdw[f] == -1:
					glColor3d(0.0,0.0,0.0)
				elif self.vdw[f] == 0:
					glColor3d(1.0,0.0,0.0)
				elif self.vdw[f] == 1:
					glColor3d(0.0,1.0,0.0)
				elif self.vdw[f] == 2:
					glColor3d(0.0,0.0,1.0)
				elif self.vdw[f] == 3:
					glColor3d(1.0,1.0,0.0)
				elif self.vdw[f] == 4:
					glColor3d(0.0,1.0,1.0)
				elif self.vdw[f] == 5:
					glColor3d(1.0,0.0,1.0)
				elif self.vdw[f] == 6:
					glColor3d(1.0,1.0,1.0)
				elif self.vdw[f] == 7:
					glColor3d(0.5,0.5,0.5)


				glNormal3d(norm1[0], norm1[1], norm1[2])
				glVertex3d(n1[0], n1[1], n1[2])
				glNormal3d(norm2[0], norm2[1], norm2[2])
				glVertex3d(n2[0], n2[1], n2[2])
				glNormal3d(norm3[0], norm3[1], norm3[2])
				glVertex3d(n3[0], n3[1], n3[2])
			glEnd()

			return

		if display_flags['show_material'] == 1 and self.id_num == display_flags['selected_index']:

			if self.state == "STATIC":
				return

			# Scan materials quickly to determine range of values (maybe not just shear modulus in future!)
			relevent_param = []
			fcolor = []

			max_val = 0
			min_val = 0
			for f in self.surface:
				#relevent_param.append(self.mat.element[f[0]].shear_modulus)
				relevent_param.append(self.mat.element[f[0]][3])
			max_val = max(relevent_param)
			min_val = min(relevent_param)

			for j in range(self.num_surface_faces):
				try:
					frac = (relevent_param[j] - min_val) / (max_val - min_val)
				except(ZeroDivisionError):
					frac = 1.0

				# Green colour map
				#fcolor.append([(50 / 255.0) * (2 - frac), 1.0, 50 * (3 - frac) / 255.0])
				#fcolor.append([(50 / 255.0) * (2 - frac),0,0])

				# Red to blue color map
				fcolor.append([1-frac,0,frac])

			# Now draw some triangles
			glBegin(GL_TRIANGLES)
	
			for f in self.surface:

				n1 = self.frames[i].node_list[f[1]][0:3]
				n2 = self.frames[i].node_list[f[2]][0:3]
				n3 = self.frames[i].node_list[f[3]][0:3]

				norm1 = self.frames[i].normal_list[f[1]]
				norm2 = self.frames[i].normal_list[f[2]]
				norm3 = self.frames[i].normal_list[f[3]]

				findex = self.surface.index(f)

				glColor3d(fcolor[findex][0], fcolor[findex][1], fcolor[findex][2])
				glNormal3d(norm1[0], norm1[1], norm1[2])
				glVertex3d(n1[0], n1[1], n1[2])
				glNormal3d(norm2[0], norm2[1], norm2[2])
				glVertex3d(n2[0], n2[1], n2[2])
				glNormal3d(norm3[0], norm3[1], norm3[2])
				glVertex3d(n3[0], n3[1], n3[2])
			glEnd()
			return

		if display_flags['binding_site_edit_mode'] == 1 and self.id_num == display_flags['selected_index']:
			glBegin(GL_TRIANGLES)

			faces_dealt_with = []
			for j in range(self.num_binding_sites):
	
				# Get type
				site_type = self.binding_site_type[j][0]

				# All faces on site
				for k in range(len(self.binding_site[j])):
					
					f = self.binding_site[j][k]
					faces_dealt_with.append(f)

					n1 = self.frames[i].node_list[self.surface[f][1]][0:3]
					n2 = self.frames[i].node_list[self.surface[f][2]][0:3]
					n3 = self.frames[i].node_list[self.surface[f][3]][0:3]

					norm1 = self.frames[i].normal_list[self.surface[f][1]]
					norm2 = self.frames[i].normal_list[self.surface[f][2]]
					norm3 = self.frames[i].normal_list[self.surface[f][3]]

					if site_type == 0:
						glColor3d(1.0,0.0,0.0)
					elif site_type == 1:
						glColor3d(0.0,1.0,0.0)
					elif site_type == 2:
						glColor3d(0.0,0.0,1.0)
					elif site_type == 3:
						glColor3d(1.0,1.0,0.0)
					elif site_type == 4:
						glColor3d(0.0,1.0,1.0)
					elif site_type == 5:
						glColor3d(1.0,0.0,1.0)
					elif site_type == 6:
						glColor3d(1.0,1.0,1.0)
					elif site_type == 7:
						glColor3d(0.5,0.5,0.5)


					glNormal3d(norm1[0], norm1[1], norm1[2])
					glVertex3d(n1[0], n1[1], n1[2])
					glNormal3d(norm2[0], norm2[1], norm2[2])
					glVertex3d(n2[0], n2[1], n2[2])
					glNormal3d(norm3[0], norm3[1], norm3[2])
					glVertex3d(n3[0], n3[1], n3[2])
			
			# Now remainder of faces		
			for f in range(self.num_surface_faces):

				if f in faces_dealt_with:
					continue
				else:
	
					n1 = self.frames[i].node_list[self.surface[f][1]][0:3]
					n2 = self.frames[i].node_list[self.surface[f][2]][0:3]
					n3 = self.frames[i].node_list[self.surface[f][3]][0:3]

					norm1 = self.frames[i].normal_list[self.surface[f][1]]
					norm2 = self.frames[i].normal_list[self.surface[f][2]]
					norm3 = self.frames[i].normal_list[self.surface[f][3]]
					
					glColor3d(0.0, 0.0, 0.0)
					glNormal3d(norm1[0], norm1[1], norm1[2])
					glVertex3d(n1[0], n1[1], n1[2])
					glNormal3d(norm2[0], norm2[1], norm2[2])
					glVertex3d(n2[0], n2[1], n2[2])
					glNormal3d(norm3[0], norm3[1], norm3[2])
					glVertex3d(n3[0], n3[1], n3[2])
			glEnd()

			return

		if display_flags['show_vdw_only'] == 1:
			glBegin(GL_TRIANGLES)
			for f in range(self.num_surface_faces):
				if self.vdw[f] == -1:
					continue
				else:
					glColor3d(1.0,1.0,0.5)

				n1 = self.frames[i].node_list[self.surface[f][1]][0:3]
				n2 = self.frames[i].node_list[self.surface[f][2]][0:3]
				n3 = self.frames[i].node_list[self.surface[f][3]][0:3]

				norm1 = self.frames[i].normal_list[self.surface[f][1]]
				norm2 = self.frames[i].normal_list[self.surface[f][2]]
				norm3 = self.frames[i].normal_list[self.surface[f][3]]

				glNormal3d(norm1[0], norm1[1], norm1[2])
				glVertex3d(n1[0], n1[1], n1[2])
				glNormal3d(norm2[0], norm2[1], norm2[2])
				glVertex3d(n2[0], n2[1], n2[2])
				glNormal3d(norm3[0], norm3[1], norm3[2])
				glVertex3d(n3[0], n3[1], n3[2])
			glEnd()

			return

		bc = display_flags['blob_colour']
		if self.blob_index == display_flags['selected_blob']:
			bc = [0,0,204/255.0]
		else:
			bc = self.normalcolor

		if display_flags['show_solid'] == 1:
			glBegin(GL_TRIANGLES)
			glColor3d(bc[0], bc[1], bc[2])
			for f in range(self.num_surface_faces):
				if self.hidden_face[f] == 1:
					continue

				n1 = self.frames[i].node_list[self.surface[f][1]][0:3]
				n2 = self.frames[i].node_list[self.surface[f][2]][0:3]
				n3 = self.frames[i].node_list[self.surface[f][3]][0:3]
				#if self.blob_index == 0:
				#	print "Frame ", i, " Blob ", self.blob_index, " Conformation ", self.conformation_index
				#	print n1, n2, n3
				#if f == 0:
				#	print n1
				#	print display_flags['selected_blob'], display_flags['selected_conformation']
				norm1 = self.frames[i].normal_list[self.surface[f][1]]
				norm2 = self.frames[i].normal_list[self.surface[f][2]]
				norm3 = self.frames[i].normal_list[self.surface[f][3]]

				glNormal3d(norm1[0], norm1[1], norm1[2])
				glVertex3d(n1[0], n1[1], n1[2])
				glNormal3d(norm2[0], norm2[1], norm2[2])
				glVertex3d(n2[0], n2[1], n2[2])
				glNormal3d(norm3[0], norm3[1], norm3[2])
				glVertex3d(n3[0], n3[1], n3[2])
			glEnd()

		if display_flags['show_mesh_surf'] == 1:
			for f in xrange(self.num_surface_faces):
				n1a = self.frames[i].node_list[self.surface[f][1]]
				n2a = self.frames[i].node_list[self.surface[f][2]]
				n3a = self.frames[i].node_list[self.surface[f][3]]

				n1 = n1a[0:3]
				n2 = n2a[0:3]
				n3 = n3a[0:3]

				norm1 = self.frames[i].normal_list[self.surface[f][1]]
				norm2 = self.frames[i].normal_list[self.surface[f][2]]
				norm3 = self.frames[i].normal_list[self.surface[f][3]]

				glBegin(GL_LINE_STRIP)	
				glColor3d(0.3,0.3,1.0)
				glNormal3d(norm1[0], norm1[1], norm1[2])
				glVertex3d(n1[0], n1[1], n1[2])
				glNormal3d(norm2[0], norm2[1], norm2[2])
				glVertex3d(n2[0], n2[1], n2[2])
				glNormal3d(norm3[0], norm3[1], norm3[2])
				glVertex3d(n3[0], n3[1], n3[2])
				glEnd()


		if display_flags['show_mesh'] == 1:
			if self.no_topology == False:
				if self.do_Fij == False:
					for el in xrange(self.num_elements):
						# Get the indices of the 4 nodes of this tetrahedron
						i1 = self.topology[el][0]
						i2 = self.topology[el][1]
						i3 = self.topology[el][2]
						i4 = self.topology[el][3]
	
						# Get the nodes
						n1a = self.frames[i].node_list[i1]
						n2a = self.frames[i].node_list[i2]
						n3a = self.frames[i].node_list[i3]
						n4a = self.frames[i].node_list[i4]
						n1 = n1a[0:3]
						n2 = n2a[0:3]
						n3 = n3a[0:3]
						n4 = n4a[0:3]
		
						# Get the surface normals at each node
						norm1 = self.frames[i].normal_list[i1]
						norm2 = self.frames[i].normal_list[i2]
						norm3 = self.frames[i].normal_list[i3]
						norm4 = self.frames[i].normal_list[i4]
	
						glBegin(GL_LINE_STRIP);
						glColor3d(0.0,0.0,0.0)
						glNormal3d(norm1[0], norm1[1], norm1[2])
						glVertex3d(n1[0], n1[1], n1[2]);
						glNormal3d(norm2[0], norm2[1], norm2[2])
						glVertex3d(n2[0], n2[1], n2[2]);
						glNormal3d(norm3[0], norm3[1], norm3[2])
						glVertex3d(n3[0], n3[1], n3[2]);
						glNormal3d(norm1[0], norm1[1], norm1[2])
						glVertex3d(n1[0], n1[1], n1[2]);
						glNormal3d(norm4[0], norm4[1], norm4[2])
						glVertex3d(n4[0], n4[1], n4[2]);
						glNormal3d(norm3[0], norm3[1], norm3[2])
						glVertex3d(n3[0], n3[1], n3[2]);
						glNormal3d(norm1[0], norm1[1], norm1[2])
						glVertex3d(n1[0], n1[1], n1[2]);
						glEnd();
				
						glBegin(GL_LINES);
						glNormal3d(norm2[0], norm2[1], norm2[2])
						glVertex3d(n2[0], n2[1], n2[2]);
						glNormal3d(norm4[0], norm4[1], norm4[2])
						glVertex3d(n4[0], n4[1], n4[2]);
						glEnd();
				else:
					for el in xrange(self.num_elements):
						# Get the indices of the 4 nodes of this tetrahedron
						i1 = self.topology[el][0]
						i2 = self.topology[el][1]
						i3 = self.topology[el][2]
						i4 = self.topology[el][3]
	
						# Get the nodes
						n1a = self.frames[i].node_list[i1]
						n2a = self.frames[i].node_list[i2]
						n3a = self.frames[i].node_list[i3]
						n4a = self.frames[i].node_list[i4]
						n1 = n1a[0:3]
						n2 = n2a[0:3]
						n3 = n3a[0:3]
						n4 = n4a[0:3]
		
						# Get the surface normals at each node
						norm1 = self.frames[i].normal_list[i1]
						norm2 = self.frames[i].normal_list[i2]
						norm3 = self.frames[i].normal_list[i3]
						norm4 = self.frames[i].normal_list[i4]
	
						# get J
						J = self.get_J(n1, n2, n3, n4)

						# get Fij and related bs
						J_inv_0 = self.first_frame_J_inv[el]
						Fij = self.mat_mult(J, J_inv_0)
						detF = self.get_det(Fij)
						FijFij = self.get_double_contraction(Fij, Fij)
						
						# get strain energy
						G = 120.0e6
						K = 640.0e6
						B = K + G/3.0
						alpha = 1.0 + G/B
						W = 1.0/(2.0 * detF) * (G * FijFij + B * (detF - alpha) * (detF - alpha) - 3 * G - B * (G/B) * (G/B))


#						c = W / self.energy_thresh


#						glBegin(GL_LINE_STRIP);
#						glColor3d(1.0,1.0,0.0)
#						glNormal3d(norm1[0], norm1[1], norm1[2])
#						glVertex3d(n1[0], n1[1], n1[2]);
#						glNormal3d(norm2[0], norm2[1], norm2[2])
#						glVertex3d(n2[0], n2[1], n2[2]);
#						glNormal3d(norm3[0], norm3[1], norm3[2])
#						glVertex3d(n3[0], n3[1], n3[2]);
#						glNormal3d(norm1[0], norm1[1], norm1[2])
#						glVertex3d(n1[0], n1[1], n1[2]);
#						glNormal3d(norm4[0], norm4[1], norm4[2])
#						glVertex3d(n4[0], n4[1], n4[2]);
#						glNormal3d(norm3[0], norm3[1], norm3[2])
#						glVertex3d(n3[0], n3[1], n3[2]);
#						glNormal3d(norm1[0], norm1[1], norm1[2])
#						glVertex3d(n1[0], n1[1], n1[2]);
#						glEnd();
#				
#						glBegin(GL_LINES);
#						glNormal3d(norm2[0], norm2[1], norm2[2])
#						glVertex3d(n2[0], n2[1], n2[2]);
#						glNormal3d(norm4[0], norm4[1], norm4[2])
#						glVertex3d(n4[0], n4[1], n4[2]);
#						glEnd();

						if W > self.energy_thresh:
							c = 1.0
						else:
							continue


						if display_flags['show_solid'] == 0:
							glDisable(GL_CULL_FACE)
							glDisable(GL_LIGHTING)

							glBegin(GL_TRIANGLES);
							glColor3d(c,0.3,0.3)

							glNormal3d(norm1[0], norm1[1], norm1[2])
							glVertex3d(n1[0], n1[1], n1[2]);
							glNormal3d(norm2[0], norm2[1], norm2[2])
							glVertex3d(n2[0], n2[1], n2[2]);
							glNormal3d(norm3[0], norm3[1], norm3[2])
							glVertex3d(n3[0], n3[1], n3[2]);

							glNormal3d(norm1[0], norm1[1], norm1[2])
							glVertex3d(n1[0], n1[1], n1[2]);
							glNormal3d(norm3[0], norm3[1], norm3[2])
							glVertex3d(n3[0], n3[1], n3[2]);
							glNormal3d(norm4[0], norm4[1], norm4[2])
							glVertex3d(n4[0], n4[1], n4[2]);

							glNormal3d(norm3[0], norm3[1], norm3[2])
							glVertex3d(n3[0], n3[1], n3[2]);
							glNormal3d(norm2[0], norm2[1], norm2[2])
							glVertex3d(n2[0], n2[1], n2[2]);
							glNormal3d(norm4[0], norm4[1], norm4[2])
							glVertex3d(n4[0], n4[1], n4[2]);

							glNormal3d(norm1[0], norm1[1], norm1[2])
							glVertex3d(n1[0], n1[1], n1[2]);
							glNormal3d(norm2[0], norm2[1], norm2[2])
							glVertex3d(n2[0], n2[1], n2[2]);
							glNormal3d(norm4[0], norm4[1], norm4[2])
							glVertex3d(n4[0], n4[1], n4[2]);

							glEnd();
							glEnable(GL_CULL_FACE)
							glEnable(GL_LIGHTING)

			else:	
				for f in xrange(self.num_surface_faces):
					n1a = self.frames[i].node_list[self.surface[f][1]]
					n2a = self.frames[i].node_list[self.surface[f][2]]
					n3a = self.frames[i].node_list[self.surface[f][3]]

					n1 = n1a[0:3]
					n2 = n2a[0:3]
					n3 = n3a[0:3]

					norm1 = self.frames[i].normal_list[self.surface[f][1]]
					norm2 = self.frames[i].normal_list[self.surface[f][2]]
					norm3 = self.frames[i].normal_list[self.surface[f][3]]

					glBegin(GL_LINE_STRIP)	
					glColor3d(1.0,1.0,0.0)
					glNormal3d(norm1[0], norm1[1], norm1[2])
					glVertex3d(n1[0], n1[1], n1[2])
					glNormal3d(norm2[0], norm2[1], norm2[2])
					glVertex3d(n2[0], n2[1], n2[2])
					glNormal3d(norm3[0], norm3[1], norm3[2])
					glVertex3d(n3[0], n3[1], n3[2])
					glEnd()

		if display_flags['show_flat'] == 1:
			glDisable(GL_LIGHTING);
			glDisable(GL_FOG);
			glShadeModel(GL_FLAT)

			glBegin(GL_TRIANGLES)
			glColor3d(bc[0], bc[1], bc[2])
			for f in xrange(self.num_surface_faces):
				n1a = self.frames[i].node_list[self.surface[f][1]]
				n2a = self.frames[i].node_list[self.surface[f][2]]
				n3a = self.frames[i].node_list[self.surface[f][3]]

				n1 = n1a[0:3]
				n2 = n2a[0:3]
				n3 = n3a[0:3]

				norm1 = self.frames[i].normal_list[self.surface[f][1]]
				norm2 = self.frames[i].normal_list[self.surface[f][2]]
				norm3 = self.frames[i].normal_list[self.surface[f][3]]

				glNormal3d(norm1[0], norm1[1], norm1[2])
				glVertex3d(n1[0], n1[1], n1[2])
				glNormal3d(norm2[0], norm2[1], norm2[2])
				glVertex3d(n2[0], n2[1], n2[2])
				glNormal3d(norm3[0], norm3[1], norm3[2])
				glVertex3d(n3[0], n3[1], n3[2])
			glEnd()

			glEnable(GL_LIGHTING);
			glEnable(GL_FOG);
			glShadeModel(GL_SMOOTH)

		if display_flags['show_node_numbers'] == 1:
			glFogfv(GL_FOG_COLOR, [1.0, 0.0, 0.0])
			glFogf(GL_FOG_DENSITY, 0.02)
			if display_flags['show_linear_nodes_only'] == 0:
				for n in xrange(self.num_nodes):
					nn = (self.frames[i].node_list[n])[0:3]
					glColor3f(1.0, 1.0, 1.0)
					glRasterPos3f(nn[0], nn[1], nn[2])
					glutBitmapString(GLUT_BITMAP_HELVETICA_18, str(n));
			else:
				for n in xrange(len(self.linear_nodes_only)):
					nn = (self.frames[i].node_list[self.linear_nodes_only[n]])[0:3]
					glColor3f(1.0, 1.0, 1.0)
					glRasterPos3f(nn[0], nn[1], nn[2])
					glutBitmapString(GLUT_BITMAP_HELVETICA_18, str(self.linear_nodes_only[n]));
			glFogfv(GL_FOG_COLOR, [.1, .1, .3])
			glFogf(GL_FOG_DENSITY, 0.002)

		if display_flags['show_pinned_nodes'] == 1:
			for n in self.pinned_nodes:
				nn = (self.frames[i].node_list[n])[0:3]
				glPointSize(10.0);
				glBegin(GL_POINTS)
				glColor3f(1.0, 0.0, 0.0)
				glVertex3f(nn[0], nn[1], nn[2])
				glEnd()

		if display_flags['show_shortest_edge'] == 1:
			if self.min_length == None:
				self.find_shortest_edge(i)
			n1 = (self.frames[i].node_list[self.shortest_edge_n1])[0:3]
			n2 = (self.frames[i].node_list[self.shortest_edge_n2])[0:3]
			glBegin(GL_LINES)
			glColor3f(1.0, 0.0, 0.0)
			glVertex3f(n1[0], n1[1], n1[2])
			glVertex3f(n2[0], n2[1], n2[2])
			glEnd()
			
			ave = [.5 * (n1[0] + n2[0]), .5 * (n1[1] + n2[1]), .5 * (n1[2] + n2[2])]
			glRasterPos3f(ave[0], ave[1], ave[2])
			glutBitmapString(GLUT_BITMAP_HELVETICA_18, str(self.min_length));

		if display_flags['show_inverted'] == 1:
			if self.no_topology == False:
				inv_str = ""
				involved_nodes = []
				for el in xrange(self.num_elements):
					# Get the indices of the 4 nodes of this tetrahedron
					i1 = self.topology[el][0]
					i2 = self.topology[el][1]
					i3 = self.topology[el][2]
					i4 = self.topology[el][3]

					# Get the nodes
					n1a = self.frames[i].node_list[i1]
					n2a = self.frames[i].node_list[i2]
					n3a = self.frames[i].node_list[i3]
					n4a = self.frames[i].node_list[i4]
					n1 = n1a[0:3]
					n2 = n2a[0:3]
					n3 = n3a[0:3]
					n4 = n4a[0:3]
	
					# Get the surface normals at each node
					norm1 = self.frames[i].normal_list[i1]
					norm2 = self.frames[i].normal_list[i2]
					norm3 = self.frames[i].normal_list[i3]
					norm4 = self.frames[i].normal_list[i4]

					# get the volume of the elements
					vol = self.get_element_volume(n1, n2, n3, n4)
					if vol * self.first_frame_vol_list[el] < 0:
						inv_str += str(el) + ","
						involved_nodes += [i1, i2, i3, i4]
						glBegin(GL_LINE_STRIP);
						glColor3d(1.0,1.0,0.0)
						glNormal3d(norm1[0], norm1[1], norm1[2])
						glVertex3d(n1[0], n1[1], n1[2]);
						glNormal3d(norm2[0], norm2[1], norm2[2])
						glVertex3d(n2[0], n2[1], n2[2]);
						glNormal3d(norm3[0], norm3[1], norm3[2])
						glVertex3d(n3[0], n3[1], n3[2]);
						glNormal3d(norm1[0], norm1[1], norm1[2])
						glVertex3d(n1[0], n1[1], n1[2]);
						glNormal3d(norm4[0], norm4[1], norm4[2])
						glVertex3d(n4[0], n4[1], n4[2]);
						glNormal3d(norm3[0], norm3[1], norm3[2])
						glVertex3d(n3[0], n3[1], n3[2]);
						glNormal3d(norm1[0], norm1[1], norm1[2])
						glVertex3d(n1[0], n1[1], n1[2]);
						glEnd();
			
						glBegin(GL_LINES);
						glNormal3d(norm2[0], norm2[1], norm2[2])
						glVertex3d(n2[0], n2[1], n2[2]);
						glNormal3d(norm4[0], norm4[1], norm4[2])
						glVertex3d(n4[0], n4[1], n4[2]);
						glEnd();
				involved_nodes = list(set(involved_nodes))
				involved_nodes = ",".join(map(str, involved_nodes))
				print "Inverted elements:"
				print inv_str
				print "Involved nodes:"
				print involved_nodes
				print "---"


	def draw_pick_frame(self, i):
		if self.num_frames == 0:
			return

		if i < 0:
			i == 0
		elif i >= self.num_frames:
			i = self.num_frames - 1

		glBegin(GL_TRIANGLES)
		pick_r = 1
		pick_g = 0
		pick_b = 0
		for f in range(self.num_surface_faces):
			n1a = self.frames[i].node_list[self.surface[f][1]]
			n2a = self.frames[i].node_list[self.surface[f][2]]
			n3a = self.frames[i].node_list[self.surface[f][3]]

			n1 = n1a[0:3]
			n2 = n2a[0:3]
			n3 = n3a[0:3]

			glColor3f(pick_r/255.0, pick_g/255.0, pick_b/255.0)
			glVertex3d(n1[0], n1[1], n1[2])
			glVertex3d(n2[0], n2[1], n2[2])
			glVertex3d(n3[0], n3[1], n3[2])

			pick_r += 1
			if pick_r == 255:
				pick_r = 0
				pick_g += 1
				if pick_g == 255:
					pick_g = 0
					pick_b += 1
		glEnd()

	def set_vdw_face(self, face_index, vdw_type):
		if face_index < 0 or face_index > self.num_surface_faces:
			print "No face picked."
			return

		if self.vdw[face_index] == vdw_type:
			self.vdw[face_index] = -1
		else:
			self.vdw[face_index] = vdw_type

	def incr_vdw_face(self, face_index):
		if face_index < 0 or face_index > self.num_surface_faces:
			print "No face picked."
			return

		self.vdw[face_index] += 1
		if self.vdw[face_index] >= 8:
			self.vdw[face_index] = -1
		print "Set face", face_index, "to", self.vdw[face_index]

	def add_face_to_binding_site(self, face_index):
		if face_index < 0 or face_index > self.num_surface_faces:
			print "No face picked."
			return
		
		# New site or not?
		if self.active_binding_site == -1:
			
			# Need a new binding site
			self.binding_site.append([face_index])
			self.num_binding_sites += 1

			# Get a type
			site_type = 3 # For now
			print "New site type = " + str(site_type)
			self.binding_site_type[face_index][0] = site_type
			self.binding_site_type[face_index][1] = len(self.binding_site) - 1
			self.active_binding_site = len(self.binding_site) - 1

		else:
			if self.binding_site_type[face_index][1] != -1:

				# Selecting a new site to append to
				self.active_binding_site = self.binding_site_type[face_index][1]
				print "New site_type = " + str(self.active_binding_site)
			else:
				self.binding_site[self.active_binding_site].append(face_index)
				#self.binding_site_type[face_index][0] = site_type
				#self.binding_site_type[face_index][1] = len(self.binding_site) - 1

	def begin_new_binding_site(self, site_type, face_index):
		
		self.binding_site_type[face_index][0] = site_type
		self.binding_site_type[face_index][1] = len(self.binding_site)
		self.binding_site.append(face_index)

	def get_dimensions(self):

		dims = [[float("inf"), -1* float("inf")] for i in range(3)]

		for n in self.frames[0].node_list:
			for i in range(3):
				if n[i] < dims[i][0]:
					dims[i][0] = n[i]
				if n[i] > dims[i][1]:
					dims[i][1] = n[i]
		return dims
	
	def get_state(self):
		return self.state

	def get_frame_state(self, i):
		if self.num_frames == 0:
			return

		if i < 0 or i >= self.num_frames:
			raise "Frame out of range"

		return self.frames[i].blob_state

	def show(self):
		self.hide_blob = False

	def hide(self):
		self.hide_blob = True

	def find_shortest_edge(self, frame_i):
		# get min length
		self.min_length = float("Inf")
		self.shortest_edge_n1 = None
		self.shortest_edge_n2 = None

		for el in xrange(self.num_elements):
			# Get the indices of the 4 nodes of this tetrahedron
			i1 = self.topology[el][0]
			i2 = self.topology[el][1]
			i3 = self.topology[el][2]
			i4 = self.topology[el][3]

			# Get the nodes
			nodes = [(self.frames[frame_i].node_list[i1])[0:3],(self.frames[frame_i].node_list[i2])[0:3],(self.frames[frame_i].node_list[i3])[0:3],(self.frames[frame_i].node_list[i4])[0:3]]

			for i in range(4):
				for j in range(i+1,4):
					dx = nodes[i][0] - nodes[j][0]
					dy = nodes[i][1] - nodes[j][1]
					dz = nodes[i][2] - nodes[j][2]
					length = math.sqrt(dx * dx + dy * dy + dz * dz)
					if length < self.min_length:
						self.min_length = length
						self.shortest_edge_n1 = i
						self.shortest_edge_n2 = j
		print "Shortest edge has length", self.min_length

	def hide_unhide_face(self, face_index):
		if face_index < 0 or face_index > self.num_surface_faces:
			print "No face picked."
			return

		self.hidden_face[face_index] *= -1
		self.vdw[face_index] = -2

	def get_element_volume(self, n0, n1, n2, n3):
		J = [	[n1[0] - n0[0], n1[1] - n0[1], n1[2] - n0[2]],
			[n2[0] - n0[0], n2[1] - n0[1], n2[2] - n0[2]],
			[n3[0] - n0[0], n3[1] - n0[1], n3[2] - n0[2]]]

		DPSI2_DX = J[2][2]*J[1][1] - J[2][1]*J[1][2];
		DPSI3_DX = J[2][1]*J[0][2] - J[2][2]*J[0][1];
		DPSI4_DX = J[1][2]*J[0][1] - J[1][1]*J[0][2];

		det = J[0][0] * DPSI2_DX + J[1][0] * DPSI3_DX + J[2][0] * DPSI4_DX;
		vol = (1.0/6.0) * det;

		return vol

	def get_det(self, F):
		a = F[2][2]*F[1][1] - F[2][1]*F[1][2];
		b = F[2][1]*F[0][2] - F[2][2]*F[0][1];
		c = F[1][2]*F[0][1] - F[1][1]*F[0][2];
		det = F[0][0] * a + F[1][0] * b + F[2][0] * c;
		return det

	def get_J_inv(self, n0, n1, n2, n3):
		J = [	[n1[0] - n0[0], n1[1] - n0[1], n1[2] - n0[2]],
			[n2[0] - n0[0], n2[1] - n0[1], n2[2] - n0[2]],
			[n3[0] - n0[0], n3[1] - n0[1], n3[2] - n0[2]]]

		J_inv = [[0,0,0], [0,0,0], [0,0,0]]

		# Construct the inverse matrix
		J_inv[0][0] = J[2][2]*J[1][1] - J[2][1]*J[1][2];
		J_inv[0][1] = J[2][1]*J[0][2] - J[2][2]*J[0][1];
		J_inv[0][2] = J[1][2]*J[0][1] - J[1][1]*J[0][2];
		J_inv[1][0] = J[2][0]*J[1][2] - J[2][2]*J[1][0];
		J_inv[1][1] = J[2][2]*J[0][0] - J[2][0]*J[0][2];
		J_inv[1][2] = J[1][0]*J[0][2] - J[1][2]*J[0][0];
		J_inv[2][0] = J[2][1]*J[1][0] - J[2][0]*J[1][1];
		J_inv[2][1] = J[2][0]*J[0][1] - J[2][1]*J[0][0];
		J_inv[2][2] = J[1][1]*J[0][0] - J[1][0]*J[0][1];
		
		# calc determinant
		det = J[0][0] * J_inv[0][0] + J[1][0] * J_inv[0][1] + J[2][0] * J_inv[0][2];
		
		# divide by determinant
		det = 1.0/det;
		J_inv[0][0]*=det; J_inv[0][1]*=det; J_inv[0][2]*=det;
		J_inv[1][0]*=det; J_inv[1][1]*=det; J_inv[1][2]*=det;
		J_inv[2][0]*=det; J_inv[2][1]*=det; J_inv[2][2]*=det;

		return J_inv

	def get_J(self, n0, n1, n2, n3):
		J = [	[n1[0] - n0[0], n1[1] - n0[1], n1[2] - n0[2]],
			[n2[0] - n0[0], n2[1] - n0[1], n2[2] - n0[2]],
			[n3[0] - n0[0], n3[1] - n0[1], n3[2] - n0[2]]]
		return J

	def get_double_contraction(self, m1, m2):
		sum = 0
		for i in range(3):
			for j in range(3):
				sum += m1[i][j] * m2[i][j]
		return sum

	def mat_mult(self, A, B):
		result = [[0,0,0], [0,0,0], [0,0,0]]
		for i in range(3):
			for j in range(3):
				for k in range(3):
					result[i][j] += A[i][k] * B[k][j]
		return result


class Frame:
	def __init__(self, blob_state, node_list, normal_list, centroid_x, centroid_y, centroid_z):
		self.blob_state = blob_state
		self.node_list = node_list
		self.normal_list = normal_list
		self.centroid_x = centroid_x
		self.centroid_y = centroid_y
		self.centroid_z = centroid_z

	def translate(self, trans):
		for n in self.node_list:
			for i in range(3):
				n[i] += trans[i]

		self.calc_centroid()

	def calc_centroid(self):
		self.centroid_x = 0.0
		self.centroid_y = 0.0
		self.centroid_z = 0.0
		for n in self.node_list:
			self.centroid_x += n[0]
			self.centroid_y += n[1]
			self.centroid_z += n[2]

		self.centroid_x *= 1.0 / len(self.node_list)
		self.centroid_y *= 1.0 / len(self.node_list)
		self.centroid_z *= 1.0 / len(self.node_list)
