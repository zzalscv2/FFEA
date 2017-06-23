# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

# from OpenGL.GL import *

# from OpenGL.GLUT import *
# from OpenGL.GLU import *
import math, os, sys
import numpy as np
import StringIO, tempfile
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
		except:
			print("\nERROR: '" + c.nodes + "' could not be loaded.")
			raise
		try:
			self.surf = FFEA_surface.FFEA_surface(c.surface)
		except:
			print("\nERROR: '" + c.surface + "' could not be loaded.")
			raise
		try:
			self.vdw = FFEA_vdw.FFEA_vdw(c.vdw)
		except:
			print("\nERROR: '" + c.vdw + "' could not be loaded.")
			raise
		
		# Only necessary for dynamic blobs
		if self.motion_state == "DYNAMIC":

			# Try to load
			try:
				self.top = FFEA_topology.FFEA_topology(c.topology)
			except:
				print("\nERROR: '" + c.topology + "' could not be loaded.")
				raise
			try:
				self.mat = FFEA_material.FFEA_material(c.material)
			except:
				print("\nERROR: '" + c.material + "' could not be loaded.")
				raise
			try:
				self.stokes = FFEA_stokes.FFEA_stokes(c.stokes)
			except:
				print("\nERROR: '" + c.stokes + "' could not be loaded.")
				raise

			try:
				self.pin = FFEA_pin.FFEA_pin(c.pin)
			except:
				print("\nERROR: '" + c.pin + "' could not be loaded.")
				raise

		# Successfully loaded, but structurally incorrect (the value self.<obj>.empty determines whether we have a default object or not i.e. not specified in script)
		if (not self.node.valid): raise IOError('Something went wrong initialising nodes')	
		if (not self.surf.valid): raise IOError('Something went wrong initialising surface')
		if (not self.vdw.valid): raise IOError('Something went wrong initialising vdw')
		if self.motion_state == "DYNAMIC":
			if (not self.top.valid): raise IOError('Something went wrong initialising topology')
			if (not self.mat.valid): raise IOError('Something went wrong initialising material')
			if (not self.stokes.valid): raise IOError('Something went wrong initialising stokes')
			if (not self.pin.valid): raise IOError('Something went wrong initialising pinned nodes')

		#
		# If certain things are not loaded, we could have reduced functionality. For example, if node.valid == False, you're screwed but if top.valid == False, you could still draw the nodes and surface
		# if self.pin == False or self.stokes == False, basically no problem whatsoever
		#

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
			#print aframe.calc_centroid(), self.init_centroid

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

	def draw_frame(self, i, frameLabel, display_flags, scale = 1.0):

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
		dan = []
		dantxt = []	
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

					# Get params
					param = self.mat.element[:,paramval]

					# Build color bins
					Erange = max(param) - min(param)
					binwidth = Erange / (num_cols - 1)
					Ebin = [min(param) + j * binwidth for j in range(num_cols)]
					Ebin[-1] = np.ceil(Ebin[-1])

					# Cheat to catch the 0 / 0 below
					if binwidth == 0.0:
						binwidth = 1.0

					# Now, draw each face
					for f in self.surf.face:
						try:
							paramfrac = (param[f.elindex] - Ebin[0]) / Erange
						except:
							paramfrac = 0.0
	
						# Get position data for viewer to draw
						n1 = self.frames[i].pos[f.n[0]][0:3]
						n2 = self.frames[i].pos[f.n[1]][0:3]
						n3 = self.frames[i].pos[f.n[2]][0:3]
						norm = self.calc_normal_2(n1, n2, n3)

						# Get color bin
						colpair = [0,1]
						for j in range(num_cols - 1):
							if param[f.elindex] >= Ebin[j] and param[f.elindex] <= Ebin[j + 1]:
								colpair = [j, j + 1]
								break

						# Get fraction through color bin (will raise error for homogeneity only)
						try:
							colfrac = (param[f.elindex] - Ebin[colpair[0]]) / binwidth

						except(ZeroDivisionError):
							colfrac = 0.0

						# What color then?
						col = (colgrad[colpair[1]] - colgrad[colpair[0]]) * colfrac + colgrad[colpair[0]]
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



			sol.extend([END])
			cmd.load_cgo(sol, display_flags['system_name'] + "_" + str(self.idnum) + "_solid_" + str(self.num_loads), frameLabel)

		#
		#  Mesh      (doable usually. catch if there's no topology i.e. STATIC blob)
		#

		if display_flags['show_mesh'] != "No Mesh":

			mes.extend( [BEGIN, LINES] )
			#mes.extend([COLOR, 1.0, 1.0, 1.0])
			#mes.extend([COLOR, 0.0, 0.0, 1.0])

			# If surface mesh, draw lines for surface only, else for entire element structure
			if display_flags['show_mesh'] == "Whole Mesh" and self.top != None:
			
				# Loop through elements
				for e in xrange(self.top.num_elements):
					n1 = self.frames[i].pos[self.top.element[e].n[0]]
					n2 = self.frames[i].pos[self.top.element[e].n[1]]
					n3 = self.frames[i].pos[self.top.element[e].n[2]]
					n4 = self.frames[i].pos[self.top.element[e].n[3]]

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

			elif display_flags['show_mesh'] == "Surface Mesh" or self.top == None:

				# Loop over surface
				#mes.extend([COLOR, 0.33, 0.33, 0.33])
				for f in xrange(self.surf.num_faces):
					n1 = self.frames[i].pos[self.surf.face[f].n[0]]
					n2 = self.frames[i].pos[self.surf.face[f].n[1]]
					n3 = self.frames[i].pos[self.surf.face[f].n[2]]

		                        mes.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
					mes.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
		                        
		                        mes.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
					mes.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )
		                       
		                        mes.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )
		                        mes.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )

			mes.extend([END])
			cmd.load_cgo(mes, display_flags['system_name'] + "_" + str(self.idnum) + "_mesh_" + str(self.num_loads), frameLabel)

		#
		#  Numbers       (again, can't always do elements)
		#

		if display_flags['show_numbers'] != "No Indices":

			# Only first frame
			if frameLabel == 1:
					
				axes = np.array([[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]])
				#scale = 0.1	 # * self.scale * self.global_scale	# Maybe change me in the future to some clever function to do with the global scale? Or get rid of global scale...
                         # No, the clever function should be a function of the shortest edge.
			 # Oh yeah right, I concur
				if display_flags['show_numbers'] == 'Node Indices':
					for n in range(self.node.num_nodes):
						nn = copy.copy(self.frames[i].pos[n][0:3])
						cyl_text(numtxt,plain,nn,str(n), scale, axes=axes)
				
				elif display_flags['show_numbers'] == 'Node Indices (Linear)':
					if len(self.linear_node_list) == 0: 
						if frameLabel == 1:
							print "Linear node indices cannot be loaded onto blob ", self.idnum, " as no topology was loaded."
							if  self.motion_state != "DYNAMIC":
								print "Try editing the FFEA script, and changing motion state to DYNAMIC"

					else:
						for n in self.linear_node_list:
							nn = copy.copy(self.frames[i].pos[n][0:3])
							cyl_text(numtxt,plain,nn,str(n), scale, axes = axes)

				elif display_flags['show_numbers'] == "Element Indicies":
						
					# Catch elements (but don't mislead i.e. no numbers
					if self.top == None:
						print "No topology! Can't show element numbers for Blob ", self.bindex
					else:
						for e in range(self.top.num_elements):
							en = self.top.element[e].calc_centroid(self.frames[i])
							cyl_text(numtxt, plain, en, str(e), scale, axes=axes)
						
				elif display_flags['show_numbers'] == "Face Indices":
					for f in range(self.surf.num_faces):
						fn = self.surf.face[f].calc_centroid(self.frames[i])
						cyl_text(numtxt, plain, fn, str(f), scale, axes=axes)
				
				# Only create object if something exists to draw (some of these above routines do nothing
				if len(numtxt) != 0:
					cmd.load_cgo(numtxt, display_flags['system_name'] + "_" + str(self.idnum) + "_numbers_" + str(self.num_loads), frameLabel)               


		#
		#  Pinned Nodes (if no pinned nodes file specified, 'break')
		#

		if display_flags['show_pinned'] == 1 and self.pin != None and self.pin.num_pinned_nodes != 0:
			pin_name = display_flags['system_name'] + "_" + str(self.idnum) + "_pinned_" + str(self.num_loads)
			psa_name = "CA"
			psa_b = 20
			text = ""
			for n in self.pin.index:
				pos = copy.copy(self.frames[i].pos[n])
				text += ("ATOM %6i %4s %3s %1s%4i    %8.3f%8.3f%8.3f\n" % (n, psa_name, "FEA", "A", n, pos[0], pos[1], pos[2]))

			# load, if it's not an empty object:
			if len(text) != 0:
				cmd.read_pdbstr(text, pin_name, frameLabel)
				cmd.show("spheres", pin_name)
				cmd.color("red", pin_name)
				
		#
		# Danger Elements! Elements that will probably invert because they have <5A lengths in them. Only draw on first frame (takes ages)
		#
		if frameLabel == 1 and display_flags['show_danger'] == 1 and self.top != None:
			#danger_name = pin_name = display_flags['system_name'] + "_" + str(self.idnum) + "_danger_" + str(self.num_loads)

			# Calculate the element lengthscales and draw all < 5A
			eindex = 0
			dindex = []
			for e in self.top.element:
				if e.get_smallest_lengthscale(self.frames[i]) / self.global_scale < 5e-10:
					dindex.append(eindex)
				eindex += 1

			# Draw the mesh
			dan.extend( [BEGIN, LINES] )
			dan.extend([COLOR, 1.0, 0.0, 0.0])
			for e in dindex:
				n1 = self.frames[i].pos[self.top.element[e].n[0]]
				n2 = self.frames[i].pos[self.top.element[e].n[1]]
				n3 = self.frames[i].pos[self.top.element[e].n[2]]
				n4 = self.frames[i].pos[self.top.element[e].n[3]]
				
				dan.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
				dan.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )

				dan.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
				dan.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )

				dan.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )
				dan.extend( [ VERTEX, n4[0], n4[1], n4[2] ] )

				dan.extend( [ VERTEX, n4[0], n4[1], n4[2] ] )
	                        dan.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )

	                        dan.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
				dan.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )

				dan.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
	                        dan.extend( [ VERTEX, n4[0], n4[1], n4[2] ] )

			dan.extend([END])
			if len(dan) != 7:
				cmd.load_cgo(dan, display_flags['system_name'] + "_" + str(self.idnum) + "_danger_" + str(self.num_loads), frameLabel)

			axes = np.array([[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]])
			
			# And the indices
			for e in dindex:
				en = self.top.element[e].calc_centroid(self.frames[i])
				cyl_text(dantxt, plain, en, str(e), scale, axes=axes)

			if len(dantxt) != 0:
				cmd.load_cgo(dantxt, display_flags['system_name'] + "_" + str(self.idnum) + "_dangernum_" + str(self.num_loads), frameLabel)

		#
		#  Load SFA: Supportive Fake Atoms #
		#		# uses read_pdbstr instead of pseudoatoms, as it is much faster!
		#
		if display_flags['load_sfa'] != "None":
			sfa_name = display_flags['system_name'] + "_" + str(self.idnum)
			psa_name = "CA"
			psa_b = 20
			text = ""
			if display_flags['load_sfa'] == "Onto Nodes":
				sfa_name += "_nfa"
				if (self.node.num_nodes > 10000):
					print "Use (ATOM number -1) to identify the node"
				for n in range(self.node.num_nodes):
					pos = (self.frames[i].pos[n].tolist())[0:3]
					if n < 10000: 
						text += ("ATOM %6i %4s %3s %1s%4i    %8.3f%8.3f%8.3f\n" % (n, psa_name, "FEA", "A", n, pos[0], pos[1], pos[2]))
					else:
						text += ("ATOM %6i %4s %3s %1s%4i    %8.3f%8.3f%8.3f\n" % (n, psa_name, "FEA", "A", 9999, pos[0], pos[1], pos[2]))


			elif display_flags['load_sfa'] == 'Onto Linear Nodes':
				if len(self.linear_node_list) == 0: 
					if frameLabel == 1:
						print "Atoms cannot be loaded onto blob ", self.idnum, " as no topology was loaded."
						if  self.motion_state != "DYNAMIC":
							print "Try editing the FFEA script, and changing motion state to DYNAMIC"

				else:
					sfa_name += "_lnfa"
					for n in self.linear_node_list:
						nn = (self.frames[i].pos[n])[0:3]
						if n == 10000: 
							print "Cannot load more than 10000 Supportive Fake Atoms"
							break
						text += ("ATOM %6i %4s %3s %1s%4i    %8.3f%8.3f%8.3f\n" % (n, psa_name, "FEA", "A", n, nn[0], nn[1], nn[2]))
	

			elif display_flags['load_sfa'] == "Onto Faces":
				sfa_name += "_ffa"
				for f in range(self.surf.num_faces):
					fn = self.surf.face[f].calc_centroid(self.frames[i])
					if f == 10000: 
						print "Cannot load more than 10000 Supportive Fake Atoms"
						break
					text += ("ATOM %6i %4s %3s %1s%4i    %8.3f%8.3f%8.3f\n" % (f, psa_name, "FEA", "A", f, fn[0], fn[1], fn[2]))


			
			elif display_flags['load_sfa'] == "Onto Elements":
				sfa_name += "_efa"
				# Catch elements (but don't mislead i.e. no numbers
				if self.top == None:
					print "No topology! Can't add atoms on elements for Blob ", self.bindex
				else:
					for e in range(self.top.num_elements):
						en = self.top.element[e].calc_centroid(self.frames[i])
						if e == 10000: 
							print "Cannot load more than 10000 Supportive Fake Atoms"
							break
						text += ("ATOM %6i %4s %3s %1s%4i    %8.3f%8.3f%8.3f\n" % (e, psa_name, "FEA", "A", e, en[0], en[1], en[2]))
						
			# in any case:
 			if len(text) > 0: 
 				cmd.read_pdbstr(text, sfa_name, frameLabel)
				cmd.hide("everything", sfa_name)
				cmd.show("spheres", sfa_name)
	

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
