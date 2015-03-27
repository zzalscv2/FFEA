from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import math

class Blob:
	def __init__(self):
		self.num_elements = 0
		self.topology = []
		self.no_topology = False
		self.num_surface_faces = []
		self.num_nodes = 0
		self.num_surface_faces = 0
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

	def load(self, idnum, nodes_fname, top_fname, surf_fname, vdw_fname, scale, blob_state, blob_pinned):
		self.id_num = idnum
		self.nodes_fname = nodes_fname
		self.nodes_scale = scale
		self.state = blob_state
		if top_fname == None:
			print "No topology file provided."
			self.no_topology = True
		else:
			self.load_topology(top_fname)
			self.no_topology = False

		self.load_surface(surf_fname)

		if vdw_fname == None:
			print "No vdw file provided. Creating a zero array."
			self.vdw = [-1 for i in xrange(self.num_surface_faces)]
		else:
			self.load_vdw(vdw_fname)

		self.hidden_face = [-1 for i in xrange(self.num_surface_faces)]
		for i in xrange(self.num_surface_faces):
			if self.vdw[i] == -2:
				self.hidden_face[i] = 1

		if blob_pinned == None:
			print "No pinned nodes file provided."
			self.num_pinned_nodes = 0
		else:
			self.load_pinned_nodes(blob_pinned)

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
		
		# skip blob, step line
		line = traj_file.readline()
#		line = traj_file.readline()
#		blob_state = line.rstrip()
		blob_state = "DYNAMIC"
		line = traj_file.readline()
		self.num_nodes = int(line)

		scale = 1e10
		nodes = []
		centroid_x = 0.0
		centroid_y = 0.0
		centroid_z = 0.0
		cdef int n
		for n in xrange(self.num_nodes):
			line = traj_file.readline().split()
			el_nodes = [float(line[i])*scale for i in xrange(10)]
			nodes.append(el_nodes)
		
			centroid_x += el_nodes[0]
			centroid_y += el_nodes[1]
			centroid_z += el_nodes[2]
		centroid_x /= self.num_nodes
		centroid_y /= self.num_nodes
		centroid_z /= self.num_nodes
		
		# Calculate average normal at each node (for gl lighting effects)
		normal_list = [[0.0, 0.0, 0.0] for i in xrange(self.num_nodes)]
		cdef int f
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
			print "Error: nodes file " + nodes_file_fname + " missing 'surface nodes:' line"
			return

		print "Scaling by " + str(self.nodes_scale)

		nodes = []
		centroid_x = 0.0
		centroid_y = 0.0
		centroid_z = 0.0
		for n in xrange(self.num_nodes):
			line = nodes_file.readline().split()
			if "interior" in line[0]:
				print "Skipping 'interior nodes:' line"
				line = nodes_file.readline().split()
			el_nodes = [float(line[i])*self.nodes_scale for i in xrange(3)]
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

	#def calc_normal(self, n1, n2, n3):
	#	ax = n2[0] - n1[0]
	#	ay = n2[1] - n1[1]
	#	az = n2[2] - n1[2]
	#	bx = n3[0] - n1[0]
	#	by = n3[1] - n1[1]
	#	bz = n3[2] - n1[2]
#
#		return [az * by - ay * bz, ax * bz - az * bx, ay * bx - ax * by]

	def get_centroid(self, i):
		if self.num_frames == 0:
			return 0.0, 0.0, 0.0

		if i < 0:
			i == 0
		elif i >= self.num_frames:
			i = self.num_frames - 1

		f = self.frames[i]
		return f.centroid_x, f.centroid_y, f.centroid_z

	def draw_frame(self, i, display_flags, iii):
		if self.hide_blob == True:
			return

		if self.num_frames == 0:
			return

		if i < 0:
			i == 0
		elif i >= self.num_frames:
			i = self.num_frames - 1

		if display_flags['hide_frozen'] == 1:
			if self.get_frame_state(i) == "FROZEN":
				return
		
		cdef int f

		if display_flags['vdw_edit_mode'] == 1 and self.id_num == display_flags['selected_blob']:
			glBegin(GL_TRIANGLES)
			for f in range(self.num_surface_faces):

				n1 = self.frames[i].node_list[self.surface[f][1]][0:3]
				n2 = self.frames[i].node_list[self.surface[f][2]][0:3]
				n3 = self.frames[i].node_list[self.surface[f][3]][0:3]

				norm1 = self.frames[i].normal_list[self.surface[f][1]]
				norm2 = self.frames[i].normal_list[self.surface[f][2]]
				norm3 = self.frames[i].normal_list[self.surface[f][3]]

				if self.vdw[f] == -1:
					glColor3d(0.3,0.3,0.0)
				else:
					glColor3d(1.0,1.0,0.5)

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
		if self.id_num == display_flags['selected_blob']:
			bc = [0.0, 1.0, 1.0]

		if display_flags['show_solid'] == 1:
			glBegin(GL_TRIANGLES)
			glColor3d(bc[0], bc[1], bc[2])
			for f in range(self.num_surface_faces):
				if self.hidden_face[f] == 1:
					continue

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

		if display_flags['show_mesh'] == 1:
			if self.no_topology == False:
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
			glFogfv(GL_FOG_COLOR, [1.0, 1.0, 1.0])
			glFogf(GL_FOG_DENSITY, 0.02)
			if display_flags['show_linear_nodes_only'] == 0:
				for n in xrange(self.num_nodes):
					nn = (self.frames[i].node_list[n])[0:3]
					glColor3f(1.0, 0.0, 0.0)
					glRasterPos3f(nn[0], nn[1], nn[2])
					glutBitmapString(GLUT_BITMAP_HELVETICA_18, str(n));
			else:
				for n in xrange(len(self.linear_nodes_only)):
					nn = (self.frames[i].node_list[self.linear_nodes_only[n]])[0:3]
					glColor3f(1.0, 0.0, 0.0)
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

class Frame:
	def __init__(self, blob_state, node_list, normal_list, centroid_x, centroid_y, centroid_z):
		self.blob_state = blob_state
		self.node_list = node_list
		self.normal_list = normal_list
		self.centroid_x = centroid_x
		self.centroid_y = centroid_y
		self.centroid_z = centroid_z
