# from OpenGL.GL import *

# from OpenGL.GLUT import *
# from OpenGL.GLU import *
import math, os, sys
import numpy as np
import FFEA_node, FFEA_surface, FFEA_topology, FFEA_material
import FFEA_stokes, FFEA_vdw, FFEA_pin, FFEA_binding_sites
import FFEA_frame

from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain

import copy

# from pymol.callback import Callback

def get_vdw_colour(index):
	if index == -1:
		return [0.5,0.5,0.5]
	elif index == 0:
		return [0.0,1.0,0.0]
	elif index == 1:
		return [1.0,0.0,0.0]
	elif index == 2:
		return [0.0,0.0,1.0]
	elif index == 3:
		return [1.0,1.0,0.0]
	elif index == 4:
		return [0.0,1.0,1.0]
	elif index == 5:
		return [1.0,0.0,1.0]
	elif index == 6:
		return [1.0,0.0,1.0]
	elif index == 7:
		return [0.5,0.0,0.0]
class Blob:
	
	def __init__(self, energy_thresh=1.06e6):
	
		self.energy_thresh = energy_thresh
		
		# Viewer flags and ids
		self.idnum = -1
		self.bindex = -1
		self.cindex = -1
		self.hide_blob = False
		self.hidden_face = None
		self.calculated_first_frame_volumes = False
		self.calculated_first_frame_J_inv = False
		self.do_Fij = False
		self.normalcolor = [32 / 255.0, 178 / 255.0, 170 / 255.0]
		
		# Structure
		self.motion_state = "DYNAMIC"
		self.top = None
		self.linear_node_list = []
		
		self.node = None
		self.surf = None
		self.mat = None
		self.vdw = None
		self.stokes = None
		self.pin = None
		self.bsites = None
		self.init_centroid = None
		self.init_rotation = None
		self.offset = np.array([0.0,0.0,0.0])
		self.min_length = None
		self.scale = 1.0
		self.global_scale = 1.0

		self.frames = []
		self.num_frames = 0
		
		self.display_flags = None

	def load(self, idnum, bindex, cindex, script):
	
		self.idnum = idnum
		self.bindex = bindex
		self.cindex = cindex
		
		b = script.blob[bindex]
		c = b.conformation[cindex]
		
		self.motion_state = c.motion_state
		
		# All will be present

		# Try to load
		try:
			self.node = FFEA_node.FFEA_node(c.nodes)
			self.surf = FFEA_surface.FFEA_surface(c.surface)
			self.vdw = FFEA_vdw.FFEA_vdw(c.vdw)
		except:
			raise

		# Successfully loaded, but structurally incorrect
		if (not self.node.valid): raise IOError('Something went wrong initialising nodes')	
		if (not self.surf.valid): raise IOError('Something went wrong initialising surface')
		if (not self.vdw.valid): raise IOError('Something went wrong initialising vdw')
		
		# only necessary for dynamic blobs
		if self.motion_state == "DYNAMIC":

			# Try to load
			try:
				self.top = FFEA_topology.FFEA_topology(c.topology)
				self.mat = FFEA_material.FFEA_material(c.material)
				self.stokes = FFEA_stokes.FFEA_stokes(c.stokes)
				self.pin = FFEA_pin.FFEA_pin(c.pin)
			except:
				raise

			# Successfully loaded, but structurally incorrect
			if (not self.top.valid): raise IOError('Something went wrong initialising topology')
			if (not self.mat.valid): raise IOError('Something went wrong initialising material')
			if (not self.stokes.valid): raise IOError('Something went wrong initialising stokes')
			if (not self.pin.valid): raise IOError('Something went wrong initialising pinned nodes')
		
		# Only necessary if kinetics are active
		if script.params.calc_kinetics == 1 and c.bsites != "":
			self.bsites = FFEA_binding_sites.FFEA_binding_sites(c.bsites)
			if (not self.bsites.valid): raise IOError('Something went wrong initialising binding sites')
		
		#
		# Calculating linear nodes only
		#
		
		if self.top != None:
			for e in self.top.element:
				for n in e.n[0:4]:
					self.linear_node_list.append(n)
		
		else:
			# Surface file uses the secondary nodes for the interactions, so it can't be used to determine the linearity
			print "Linear nodes cannot be known without a topology."

		self.linear_node_list = set(self.linear_node_list)
		
		# Any initialisation done in ffea?
		if b.centroid != None:
			try:
				self.init_centroid = np.array(b.centroid)	
			except:
				self.init_centroid = None

		if b.rotation != None:		
			try:
				self.init_rotation = np.array(b.rotation)
			except:
				self.init_rotation = None
				
		try:
			self.scale = b.scale
		except:
			self.scale = 1.0
		
		# Initialise stuff that we didn't get
		self.hidden_face = [-1 for i in range(self.surf.num_faces)]
		if self.vdw != None and self.vdw.num_faces == 0:
			self.vdw.set_num_faces(self.surf.num_faces)
			
		# Calculate stuff that needs calculating
	
	'''	
	def load(self, idnum, blob_index, conformation_index, nodes_fname, top_fname, surf_fname, vdw_fname, scale, blob_state, blob_pinned, blob_mat, binding_fname, blob_centroid_pos, blob_rotation, ffea_fpath = "."):
		self.idnum = idnum
		self.bindex = blob_index
		self.cindex = conformation_index
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
		print self.state
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
			self.vdw = [-1 for i in xrange(self.surf.num_faces)]
		else:
                        if os.path.isabs(vdw_fname) == False:
                             vdw_fname = os.path.join(ffea_fpath, vdw_fname)
			try:
				self.load_vdw(vdw_fname)
			except(IOError):
				print "Vdw file '" + vdw_fname + "' not found. Creating a zero array."
				self.vdw = [-1 for i in xrange(self.surf.num_faces)]
			except:
				print "Unknow error in 'self.load_vdw()'. Creating a zero array."
				self.vdw = [-1 for i in xrange(self.surf.num_faces)]

		self.hidden_face = [-1 for i in xrange(self.surf.num_faces)]
		for i in xrange(self.surf.num_faces):
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
	'''
	def load_topology(self, top_fname):
		print "Reading in topology file " + top_fname
       
		top = open(top_fname, "r")
		line = top.readline().rstrip()
		if line != "ffea topology file" and line != "walrus topology file":
			print "Error: Topology file " + top_fname + " missing 'ffea topology file' first line"
			return

		line = top.readline().split()
		self.elem.num_elements = int(line[1])
		print "num_elements = ", self.elem.num_elements

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
		self.surf.num_faces = int(line[1])
		print "num_surface_faces = ", self.surf.num_faces

		line = surf.readline().rstrip()
		if line != "faces:":
			print "Error: surface file " + surf_fname + " missing 'faces:' line"
			return

		self.surface = []
		for n in xrange(self.surf.num_faces):
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

		if num_vdw_faces != self.surf.num_faces:
			print "Error. Number of faces in vdw file (" + str(num_vdw_faces) + ") does not match number of faces in surface file (" + str(self.surf.num_faces) + ")"
			return

		line = vdw_file.readline().rstrip()
		if line != "vdw params:":
			print "Error: vdw file " + vdw_fname + " missing 'vdw params:' line"
			return

		self.vdw = []
		for n in xrange(self.surf.num_faces):
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
		self.binding_site_type = [[-1, -1] for i in range(self.surf.num_faces)]
		
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
		vdw_file.write("num_faces " + str(self.surf.num_faces) + "\n")
		vdw_file.write("vdw params:\n")
		for n in xrange(self.surf.num_faces):
			vdw_file.write(str(self.vdw[n]) + "\n")
		vdw_file.close()
		print "Finished writing vdw file " + vdw_fname

#	def load_frame(self, traj_file):
#
#		# Inactive conf
#		if traj_file == None:
#			self.frames.append(None)
#			self.num_frames += 1
#			return
#
#		# There is a trajectory! Ignore offset
#		self.offset = [0.0, 0.0, 0.0]
#
#		# skip blob, step line
#		line = traj_file.readline()
#		blob_state = traj_file.readline().rstrip()
#		
#		nodes = []
#		centroid_x = 0.0
#		centroid_y = 0.0
#		centroid_z = 0.0
#		#cdef int n
#		for n in xrange(self.num_nodes):
#			line = traj_file.readline().split()
#			el_nodes = [float(line[i])* self.global_scale for i in xrange(10)]
#			nodes.append(el_nodes)
#		
#			centroid_x += el_nodes[0]
#			centroid_y += el_nodes[1]
#			centroid_z += el_nodes[2]
#
#		centroid_x /= self.num_nodes
#		centroid_y /= self.num_nodes
#		centroid_z /= self.num_nodes
#		
 #               ## ## Normals are calculated later ## 
#		# # Calculate average normal at each node (for gl lighting effects)
#		# normal_list = [[0.0, 0.0, 0.0] for i in xrange(self.num_nodes)]
#		# #cdef int f
#		# for f in xrange(self.surf.num_faces):
#		#	# get node indices of this face
#		#	i1 = self.surf.face[f].n[0]
#		#	i2 = self.surf.face[f].n[1]
#		#	i3 = self.surf.face[f].n[2]
#		
#		#	# get the normal of the face
#		#	norm = self.calc_normal(nodes[i1], nodes[i2], nodes[i3])
#		
#		#	normal_list[i1][0] += norm[0]
#		#	normal_list[i1][1] += norm[1]
#		#	normal_list[i1][2] += norm[2]
#		
#		#	normal_list[i2][0] += norm[0]
#		#	normal_list[i2][1] += norm[1]
#		#	normal_list[i2][2] += norm[2]
#		
#		#	normal_list[i3][0] += norm[0]
#		#	normal_list[i3][1] += norm[1]
#		#	normal_list[i3][2] += norm[2]
#		
#		
#		# self.frames.append(Frame(blob_state, nodes, normal_list, centroid_x, centroid_y, centroid_z))
#		self.frames.append(Frame(blob_state, nodes, NULL, centroid_x, centroid_y, centroid_z))
#		self.num_frames += 1
#
#		if self.calculated_linear_nodes == False:
#			if self.no_topology == False:
#				self.linear_node_list = []
#				for i in range(self.num_nodes):
#					for el in self.topology:
#						if i in el[0:4]:
#							self.linear_node_list.append(i)
#							break
#				print "Found", len(self.linear_node_list), "linear nodes."
#				self.calculated_linear_nodes = True	
#
#		if self.calculated_first_frame_volumes == False:
#			self.first_frame_vol_list = []
#			if self.no_topology == False:
#				print "Calculating volumes of elements in first frame..."
#
#				# Find first frame available
#				for i in range(len(self.frames)):
#					if self.frames[i] != None:
#						frame_index = i
#						break
#
#				for el in xrange(self.elem.num_elements):
#					# Get the indices of the 4 nodes of this tetrahedron
#					i1 = self.topology[el][0]
#					i2 = self.topology[el][1]
#					i3 = self.topology[el][2]
#					i4 = self.topology[el][3]
#
#					# Get the nodes
#					n1a = self.frames[frame_index].pos[i1]
#					n2a = self.frames[frame_index].pos[i2]
#					n3a = self.frames[frame_index].pos[i3]
#					n4a = self.frames[frame_index].pos[i4]
#					n1 = n1a[0:3]
#					n2 = n2a[0:3]
#					n3 = n3a[0:3]
#					n4 = n4a[0:3]
#	
#					# get the volume of the elements
#					self.first_frame_vol_list.append(self.get_element_volume(n1, n2, n3, n4))
#				print "...Done."
#				self.calculated_first_frame_volumes = True
#
#
#		if self.calculated_first_frame_J_inv == False:
#			self.first_frame_J_inv = []
#			if self.no_topology == False:
#				print "Calculating J inv of elements in first frame..."
#
#				# Find first frame available
#				for i in range(len(self.frames)):
#					if self.frames[i] != None:
#						frame_index = i
#						break
#
#				for el in xrange(self.elem.num_elements):
#					# Get the indices of the 4 nodes of this tetrahedron
#					i1 = self.topology[el][0]
#					i2 = self.topology[el][1]
#					i3 = self.topology[el][2]
#					i4 = self.topology[el][3]
#
#					# Get the nodes
#					n1a = self.frames[frame_index].pos[i1]
#					n2a = self.frames[frame_index].pos[i2]
#					n3a = self.frames[frame_index].pos[i3]
#					n4a = self.frames[frame_index].pos[i4]
#					n1 = n1a[0:3]
#					n2 = n2a[0:3]
#					n3 = n3a[0:3]
#					n4 = n4a[0:3]
#	
#					# get the J_inv of the elements
#					self.first_frame_J_inv.append(self.get_J_inv(n1, n2, n3, n4))
#				print "...Done."
#				self.calculated_first_frame_J_inv = True
 

	def delete_all_frames(self):
		self.num_frames = 0
		self.frames = []

	def set_dead_frame(self):
		self.frames.append(None)
		self.num_frames += 1
	
	def set_nodes_as_frame(self):
	
		print "Setting nodes as initial frame..."
		
		# Get a frame
		aframe = FFEA_frame.FFEA_frame()
		aframe.build_from_node(self.node)
		
		# Move and rotate it
		if self.init_centroid != None:
			print "=============================="
			print "Moving to starting position..."
			print "=============================="
            
			aframe.set_pos(self.init_centroid)
			print aframe.calc_centroid()

		if self.init_rotation != None:
			print "=============================="
			print "Rotating to starting orientation..."
			print "=============================="
			aframe.rotate(self.init_rotation)

		# Now scale
		aframe.scale(self.scale * self.global_scale)

		# Append it to the list
		self.frames.append(aframe)
		self.num_frames += 1
	
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

                # vx = az * by - ay * bz
                # vy = ax * bz - az * bx
                # vz = ay * bx - ax * by
                # n = math.sqrt(vx * vx + vy * vy + vz * vz) 
                # return ( vx/n, vy/n, vz/n )
 
	def calc_normal_2(self, n1, n2, n3):
		ax = n2[0] - n1[0]
		ay = n2[1] - n1[1]
		az = n2[2] - n1[2]
		bx = n3[0] - n2[0]
		by = n3[1] - n2[1]
		bz = n3[2] - n2[2]

		return [az * by - ay * bz, ax * bz - az * bx, ay * bx - ax * by]

                # vx = az * by - ay * bz
                # vy = ax * bz - az * bx
                # vz = ay * bx - ax * by
                # n = math.sqrt(vx * vx + vy * vy + vz * vz) 
                # return ( vx/n, vy/n, vz/n )
 
	def set_num_loads(self, i):
		self.num_loads = i

	def calc_centroid(self, i):

		if self.motion_state == "STATIC":
			i = 0

		if self.num_frames == 0:
			return None, None, None

		if self.frames[i] == None:
			return None, None, None

		#if i < 0:
		#	i = 0

		elif i >= self.num_frames:
			i = self.num_frames - 1

		f = self.frames[i]
		return f.calc_centroid()

	def draw_frame(self, i, frameLabel, display_flags):

		# Make a copy of the display flags so the user input one doesn't change!
		
		# Ideally these checks shouldn't ned to be here, but whatever
		if self.motion_state == "STATIC":
			i = 0

		if self.num_frames == 0:
			print "num_frames = 0"
			return

		if self.frames[i] == None:
			print "frame list empty"
			return

		if self.hide_blob == True:
			"blob hidden"
			return

		if i < 0:
			i = self.num_frames + i
		elif i >= self.num_frames:
			i = self.num_frames - 1
		
		sol = []
		mes = []	
		numtxt = []
		pinsphere = []

		print "loading frame ", frameLabel, " for blob ", self.idnum
		frameLabel += 1

		#
		#  Solid
		#

		MatOpt = ["Density", "Shear Viscosity", "Bulk Viscosity", "Shear Modulus", "Bulk Modulus", "VdW"]
		if display_flags['matparam'] != "No Solid":
		        sol.extend( [ BEGIN, TRIANGLES ] )

			# Can we draw material properties?
			default = False
			if MatOpt.count(display_flags['matparam']) and self.mat == None:
				if display_flags['matparam'] != "VdW":
					print "Cannnot draw material params for blob " + str(self.bindex) + ". Defaulting..."
					default = True
			
			# If solid, draw all triangles
			if default or display_flags['matparam'] == "Plain Solid":
		                # sol.extend( [ COLOR, bc[0], bc[1], bc[2] ] )
				for f in range(self.surf.num_faces):
					if self.hidden_face[f] == 1:
						continue

					n1 = self.frames[i].pos[self.surf.face[f].n[0]][0:3]
					n2 = self.frames[i].pos[self.surf.face[f].n[1]][0:3]
					n3 = self.frames[i].pos[self.surf.face[f].n[2]][0:3]
				
		                        norm = self.calc_normal_2(n1, n2, n3)

		                        sol.extend( [ NORMAL, -norm[0], -norm[1], -norm[2] ] )
		                        sol.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
		                        sol.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
		                        sol.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )

			elif MatOpt.count(display_flags['matparam']) == 1:

				# material drawing
				
				# Get param
				if display_flags['matparam'] == "Density":
					paramval = 0
				elif display_flags['matparam'] == "Shear Viscosity":
					paramval = 1
				elif display_flags['matparam'] == "Bulk Viscosity":
					paramval = 2
				elif display_flags['matparam'] == "Shear Modulus":
					paramval = 3
				elif display_flags['matparam'] == "Bulk Modulus":
					paramval = 4
				elif display_flags['matparam'] == "VdW":
					paramval = 5

				if (paramval != 5): ## that means we do proper material parameters

					# Get range of colours
					colgrad = [np.array([0.0,0.0,1.0]), np.array([0.0,1.0,0.0]), np.array([1.0,1.0,0.0]), np.array([1.0,0.0,0.0])]	# Blue green yellow red
					#colgrad = [np.array([1.0,0.0,0.0]), np.array([1.0,1.0,0.0]), np.array([0.0,1.0,0.0]), np.array([0.0,0.0,1.0])]	# Red yellow green blue
					num_cols = len(colgrad)
					param = []
	
					for e in self.mat.element:
						param.append(e[paramval])
			
					# Relative or absolute max / min here
					maxval = max(param)
					minval = min(param)
					#maxval = 100e9
					#minval = 100e7
					rangeval = maxval - minval
	
					# Now, draw each face
					for f in self.surf.face:
						
						n1 = self.frames[i].pos[f.n[0]][0:3]
						n2 = self.frames[i].pos[f.n[1]][0:3]
						n3 = self.frames[i].pos[f.n[2]][0:3]
			                        norm = self.calc_normal_2(n1, n2, n3)
	
						# Calc and add colour first
						if rangeval == 0.0:
							paramfrac = 1.0
						else:
							paramfrac = (param[f.elindex] - minval) / rangeval

						# Which interval?
						if paramfrac <= 0:
							colpair = [0, 0]
						elif paramfrac >= 1:
							colpair = [num_cols - 1, num_cols - 1]
						else:
							colpair = [int(np.floor(paramfrac * num_cols)), int(np.floor(paramfrac * num_cols)) + 1]

						# Where in interval
						col = (colgrad[colpair[1]] - colgrad[colpair[0]]) * paramfrac + colgrad[colpair[0]]
						#print col[0], col[1], col[2]
						sol.extend([COLOR, col[0], col[1], col[2]])
			                        sol.extend( [ NORMAL, -norm[0], -norm[1], -norm[2] ] )
			                        sol.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
			                        sol.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
			                        sol.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )


				else: ## in that case, plot VdW! 
					for f in range(self.surf.num_faces):
						if self.hidden_face[f] == 1:
							continue
	
						n1 = self.frames[i].pos[self.surf.face[f].n[0]][0:3]
						n2 = self.frames[i].pos[self.surf.face[f].n[1]][0:3]
						n3 = self.frames[i].pos[self.surf.face[f].n[2]][0:3]
					
			                        norm = self.calc_normal_2(n1, n2, n3)
	
						if self.vdw.index[f] != self.vdw.index[f - 1] or f == 0:
							bc = get_vdw_colour(self.vdw.index[f])
							sol.extend([ COLOR, bc[0], bc[1], bc[2] ])
	
			                        sol.extend( [ NORMAL, -norm[0], -norm[1], -norm[2] ] )
			                        sol.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
			                        sol.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
			                        sol.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )



			sol.append(END)
			cmd.load_cgo(sol, display_flags['system_name'] + "_" + str(self.idnum) + "_solid_load_" + str(self.num_loads), frameLabel)

		#
		#  Mesh      (doable usually. catch if there's no topology i.e. STATIC blob)
		#

		if display_flags['show_mesh'] != "No Mesh":
	
			# Check for topology
			if self.top == None:
				 display_flags['show_mesh'] = "Surface Mesh"

			mes.extend( [BEGIN, LINES] )
			#mes.extend([COLOR, 1.0, 1.0, 1.0])
			#mes.extend([COLOR, 0.0, 0.0, 1.0])

			# If surface mesh, draw lines for surface only, else for entire element structure
			if display_flags['show_mesh'] == "Whole Mesh":
			
				# Loop through elements
				for e in xrange(self.top.num_elements):
					n1a = self.frames[i].pos[self.top.element[e].n[0]]
					n2a = self.frames[i].pos[self.top.element[e].n[1]]
					n3a = self.frames[i].pos[self.top.element[e].n[2]]
					n4a = self.frames[i].pos[self.top.element[e].n[3]]

					n1 = n1a[0:4]
					n2 = n2a[0:4]
					n3 = n3a[0:4]
					n4 = n4a[0:4]

		                        mes.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
					mes.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )

					mes.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
					mes.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )

					mes.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )
		                        mes.extend( [ VERTEX, n4[0], n4[1], n4[2] ] )

		                        mes.extend( [ VERTEX, n4[0], n4[1], n4[2] ] )
		                        mes.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )

		                        mes.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
					mes.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )

					mes.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
		                        mes.extend( [ VERTEX, n4[0], n4[1], n4[2] ] )

			elif display_flags['show_mesh'] == "Surface Mesh":

				# Loop over surface
				for f in xrange(self.surf.num_faces):
					n1a = self.frames[i].pos[self.surf.face[f].n[0]]
					n2a = self.frames[i].pos[self.surf.face[f].n[1]]
					n3a = self.frames[i].pos[self.surf.face[f].n[2]]

					n1 = n1a[0:3]
					n2 = n2a[0:3]
					n3 = n3a[0:3]


					# norm1 = self.frames[i].normal_list[self.surf.face[f].n[0]]
					# norm2 = self.frames[i].normal_list[self.surf.face[f].n[1]]
					# norm3 = self.frames[i].normal_list[self.surf.face[f].n[2]]

		                        # mes.extend( [ LINEWIDTH, 2.0 ] )
		                        # mes.extend( [ COLOR, 0.3, 0.3, 1.0 ] )
		                        # mes.extend( [ NORMAL, norm1[0], norm1[1], norm1[2] ] )
		                        mes.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
					mes.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
		                        # mes.extend( [ NORMAL, norm2[0], norm2[1], norm2[2] ] )
		                        mes.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
					mes.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )
		                        # mes.extend( [ NORMAL, norm3[0], norm3[1], norm3[2] ] )
		                        mes.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )
		                        mes.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )

			mes.append(END)
			cmd.load_cgo(mes, display_flags['system_name'] + "_" + str(self.idnum) + "_mesh_load_" + str(self.num_loads), frameLabel)

		#
		#  Numbers       (again, can't always do elements)
		#

		if display_flags['show_numbers'] != "No Indices":

			# Only first frame
			if frameLabel == 1:
					
				axes = np.array([[15.0,0.0,0.0],[0.0,15.0,0.0],[0.0,0.0,15.0]])
				scale = 0.1	 # * self.scale * self.global_scale	# Maybe change me in the future to some clever function to do with the global scale? Or get rid of global scale...
                         # No, the clever function should be a function of the shortest edge.
	
				if display_flags['show_numbers'] == 'Node Indices':
					for n in range(self.node.num_nodes):
						nn = (self.frames[i].pos[n])[0:3]
						cyl_text(numtxt,plain,nn,str(n), scale, axes=axes * scale)
	
				elif display_flags['show_numbers'] == 'Node Indices (Linear)':
					for n in self.linear_node_list:
						nn = (self.frames[i].pos[n])[0:3]
						cyl_text(numtxt,plain,nn,str(n),0.10 * scale, axes=axes * scale)
	
				elif display_flags['show_numbers'] == "Element Indicies":
						
					# Catch elements (but don't mislead i.e. no numbers
					if self.top == None:
						print "No topology! Can't show element numbers for Blob ", self.bindex
					else:
						for e in range(self.top.num_elements):
							en = self.top.element[e].calc_centroid(self.frames[i])
							cyl_text(numtxt, plain, en, str(e), 0.10 * scale, axes=axes * scale)
						
				elif display_flags['show_numbers'] == "Face Indices":
					for f in range(self.surf.num_faces):
						fn = self.surf.face[f].calc_centroid(self.frames[i])
						cyl_text(numtxt, plain, fn, str(f), 0.10 * scale, axes=axes * scale)
				
				# Only create object if something exists to draw (some of these above routines do nothing
				if len(numtxt) != 0:
					cmd.load_cgo(numtxt, display_flags['system_name'] + "_" + str(self.idnum) + "_numbers_load_" + str(self.num_loads), frameLabel)               


		#
		#  Pinned Nodes (if no pinned nodes file specified, 'break')
		#

		if display_flags['show_pinned'] == 1 and self.pin != None and self.pin.num_pinned_nodes != 0:
			
			pinsphere.extend([COLOR, 1.0, 0.0, 0.0])
			for n in self.pin.index:
				nn = self.frames[i].pos[n]
				pinsphere.extend( [ SPHERE, nn[0], nn[1], nn[2], 5 ] )

			# Only load if there are pinned nodes (dur!)
			if len(pinsphere) != 0:
				cmd.load_cgo(pinsphere, display_flags['system_name'] + "_" + str(self.idnum) + "_pinned_load_" + str(self.num_loads), frameLabel) 


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
		for f in range(self.surf.num_faces):
			n1a = self.frames[i].pos[self.surf.face[f].n[0]]
			n2a = self.frames[i].pos[self.surf.face[f].n[1]]
			n3a = self.frames[i].pos[self.surf.face[f].n[2]]

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
		if face_index < 0 or face_index > self.surf.num_faces:
			print "No face picked."
			return

		if self.vdw[face_index] == vdw_type:
			self.vdw[face_index] = -1
		else:
			self.vdw[face_index] = vdw_type

	def incr_vdw_face(self, face_index):
		if face_index < 0 or face_index > self.surf.num_faces:
			print "No face picked."
			return

		self.vdw[face_index] += 1
		if self.vdw[face_index] >= 8:
			self.vdw[face_index] = -1
		print "Set face", face_index, "to", self.vdw[face_index]

	def add_face_to_binding_site(self, face_index):
		if face_index < 0 or face_index > self.surf.num_faces:
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

		for n in self.frames[0].pos:
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

		for el in xrange(self.elem.num_elements):
			# Get the indices of the 4 nodes of this tetrahedron
			i1 = self.topology[el][0]
			i2 = self.topology[el][1]
			i3 = self.topology[el][2]
			i4 = self.topology[el][3]

			# Get the nodes
			nodes = [(self.frames[frame_i].pos[i1])[0:3],(self.frames[frame_i].pos[i2])[0:3],(self.frames[frame_i].pos[i3])[0:3],(self.frames[frame_i].pos[i4])[0:3]]

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
		if face_index < 0 or face_index > self.surf.num_faces:
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
		try:
			det = 1.0/det;
		except(ZeroDivisionError):
			for i in range(3):
				for j in range(3):
					if i == j:
						J_inv[i][j] = 1
					else:
						J_inv[i][j] = 0
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
