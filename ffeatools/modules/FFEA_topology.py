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

import os, sys
from time import sleep
import numpy as np
import FFEA_surface
from FFEA_exceptions import *

class FFEA_topology:

	def __init__(self, fname = ""):
	
		self.reset()

		if fname == "":
			self.valid = True
			sys.stdout.write("done! Empty object initialised.\n")
			return

		try:
			self.load(fname)

		except FFEAFormatError as e:
			self.reset()
			print_error()
			print("Formatting error at line " + e.lin + "\nLine(s) should be formatted as follows:\n\n" + e.lstr)
			raise

		except FFEAIOError as e:
			self.reset()
			print_error()
			print("Input error for file " + e.fname)
			if e.fext != [""]:
				print("       Acceptable file types:")
				for ext in e.fext:
					print("       " + ext)
		except IOError:
			raise

	def load(self, fname):

		sys.stdout.write("Loading FFEA topology file...")

		# File format?
		base, ext = os.path.splitext(fname)
		try:
			if ext == ".top":
				self.load_top(fname)
			elif ext == ".ele":
				self.load_ele(fname)
			elif ext == ".vol":
				self.load_vol(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".top", ".ele", ".vol"])

		except:
			raise

		self.valid = True
		self.empty = False
		sys.stdout.write("done!\n")

	def load_top(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea topology file" and line != "walrus topology file":
			print("\tExpected 'ffea topology file' but found " + line)
			raise TypeError

		num_elements = int(fin.readline().split()[1])
		num_surface_elements = int(fin.readline().split()[1])
		num_interior_elements = int(fin.readline().split()[1])

		fin.readline()

		# Read elements now
		eltype = 0	
		while(True):
			sline = fin.readline().split()

			# Get an element
			try:
				if sline[0].strip() == "interior":
					eltype = 1
					continue
	
				elif len(sline) == 4:
					el = FFEA_element_tet_lin()
				elif len(sline) == 10:
					el = FFEA_element_tet_sec()

			except(IndexError):
				break

			el.set_indices(sline)
			self.add_element(el, eltype=eltype)

		fin.close()

	def load_vol(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Find start of elements
		while(True):
			line = fin.readline()
			if line == "":
				raise IOError
			
			elif line.strip() == "volumeelements":
				break

		# Read num_elements
		num_elements = int(fin.readline())

		# Get all elements
		zeroindexing = False
		for i in range(num_elements):
			sline = fin.readline().split()[2:]
			el = FFEA_element_tet_lin()
			
			# Test indices
			for s in sline:
				if s.strip() == "0":
					zeroindexing = True
					break
 
			el.set_indices(sline)
			self.add_element(el)

		# Indexing from 0
		if not zeroindexing:
			for i in range(self.num_elements):
				for j in range(4):
					self.element[i].n[j] -= 1

	def load_ele(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Test format
		sline = fin.readline().split()
		if len(sline) != 3:
			print("\tExpected '<num_elements> <num_vertices> 0' but found " + line)
			raise TypeError

		num_elements = int(sline[0])
		num_interior_elements = 0
		num_surface_elements = num_elements

		# Read elements now
		while(True):
			sline = fin.readline().split()
			if sline[0].strip() == "#":
				break

			sline = sline[1:]
			# Get an element
			try:
				if len(sline) == 4:
					el = FFEA_element_tet_lin()
				elif len(sline) == 10:
					el = FFEA_element_tet_sec()

			except(IndexError):
				break

			# Indexing from zero
			sline = [int(s) - 1 for s in sline]
			el.set_indices(sline)
			self.add_element(el)

		fin.close()

	def add_element(self, el, eltype = -1):

		if eltype == -1:
			self.num_surface_elements += 1
			el.interior = None

		elif eltype == 0:
			self.num_surface_elements += 1
			el.interior = False
		else:
			self.num_interior_elements += 1
			el.interior = True

		self.element.append(el)
		self.num_elements += 1

	def get_num_elements(self):
		return len(self.element)

	def get_linear_nodes(self):
	
		# Get them all
		n = []
		for e in self.element:
			for index in e.n[0:4]:
				n.append(index)
		
		# Make a list of a set
		return list(set(n))

	def calc_CoM(self, node, mat):

		CoM = np.array([0.0,0.0,0.0])
		eindex = -1
		total_mass = 0.0
		for e in self.element:
			eindex += 1
			elmass = e.calc_volume(node) * mat.element[eindex][0]
			total_mass += elmass			
			CoM += elmass * (np.mean([node.pos[n] for n in e.n[0:4]], axis=0))
		self.CoM = CoM * 1.0/total_mass
		return self.CoM

	def get_CoM(self):
		return self.CoM
	
	def extract_surface(self):
		
		faces = []
		facessorted = []
		elindex = []
		surf = FFEA_surface.FFEA_surface()

		# Get all faces in entire system. If only occurs once, it's a surface face. Keep the order
		for i in range(self.num_elements):
			e = self.element[i]

			if len(e.n) == 10:
				print("Error. Functionality not yet available for 2nd order elements")
				return

			for j in range(4):
				faces.append(e.get_linear_face(j, obj=False))
				facessorted.append(sorted(e.get_linear_face(j, obj=False)))
				elindex.append(i)

		# Now, go through sorted faces and if it occurs twice, delete both from the actual surface
		run = -1
		while(facessorted != []):
			run += 1
			# Get a face
			fsort = facessorted.pop(0)					
			f = faces.pop(0)
			eind = elindex.pop(0)

			# Is it in the list?
			if fsort not in facessorted:
				
				# It was only there once! Surface face baby
				sf = FFEA_surface.FFEA_face_tri_lin()
				sf.set_indices(f, elindex = eind)				
				surf.add_face(sf)
			else:
				# Delete the other face in the list
				for i in range(len(facessorted)):
					if facessorted[i] == fsort:
						index = i
						break
				try:
					facessorted.pop(index)
					faces.pop(index)
					elindex.pop(index)
				except:
					break
	
		return surf

	def calculateInterior(self, surf=None):

		# Don't continue without surface (we could do, but it's slow)
		if surf == None:
			print("Cannot calculate interior elements without an associated surface (currently).")
			return

		# Don't continue if we're already done
		num_done = 0
		for e in self.element:
			if e.interior == None:
				break
			else:
				num_done += 1

		if num_done == self.num_elements:
			return

		# Set all elements as default to interior elements
		self.num_interior_elements = self.num_elements
		self.num_surface_elements = 0

		for i in range(self.num_elements):
			self.element[i].interior = True


		# Now, use surface to work out which are surface elements
		for i in range(surf.num_faces):
			if self.element[surf.face[i].elindex].interior:
				self.element[surf.face[i].elindex].interior = False
				self.num_interior_elements -= 1
				self.num_surface_elements += 1

		# Get a map, so we know what element wil go where
		amap = [-1 for i in range(self.num_elements)]
		index = 0

		for i in range(self.num_elements):
			if self.element[i].interior:
				amap[i] = index 
				index += 1

		for i in range(len(amap)):
			if amap[i] == -1:
				amap[i] = index
				index += 1

		# Now, reorder actual elements
		old_els = self.element
		self.element = [None for i in range(self.num_elements)]

		for i in range(self.num_elements):
			self.element[amap[i]] = old_els[i]

		# And reorder the surface indices
		for i in range(surf.num_faces):
			surf.face[i].elindex = amap[surf.face[i].elindex]

	def isElementInterior(self, index):
		
		try:
			testEl = self.element[index]
		except(IndexError):
			print("Element ", index, " does not exist.")
			return False

		# First, see if already calculated. Else, set default assumption
		if testEl.interior != None:
			return testEl.interior

		else:
			testEl.interior = False

		# To test if interior, see if a faces are repeated. If all are, element is interior
		i = -1
		faces_not_found = [0,1,2,3]
		for el in self.element:
			
			i += 1
			el_connected = False
			faces_to_remove = -1

			# Same element doesn't count
			if i == index:
				continue
		
			for j in faces_not_found:
				testFace = testEl.get_linear_face(j)
				
				for k in range(4):
					face = el.get_linear_face(k)

					if face.isSame(testFace):
						faces_to_remove = j
						el_connected = True
						break

				if el_connected:
					break
	
			try:
				faces_not_found.remove(faces_to_remove)
			except:
				pass
			if faces_not_found == []:
				testEl.interior = True
				break

		return testEl.interior
	
	def increase_order(self, node = None, surf = None, stokes = None):
		
		# Increases the order of this topology, and all of the associated structures (if they exist)

		# Check current order (function currently only for 1st order - 2nd order)
		for e in self.element:
			if len(e.n) == 10:
				print("Error. Increasing to order > 2 is currently not supported")
				return

		# Get a unique edge list (edge is two nodes, and a potential 3rd node)
		# Also get num_nodes whilst we're at it (in case node object not provided)
		edges = []
		max_node_index = -1 * float("inf")
		for i in range(self.num_elements):
			n = self.element[i].n
			n = sorted(n)

			nmax = max(n)
			if nmax > max_node_index:
				max_node_index = nmax

			for j in range(3):
				for k in range(j + 1, 4):
					ed = (n[j], n[k], -1)
					edges.append(ed)
			
		edges = list(set(edges))

		# For each edge, add a new midpoint
		new_edges = []
		for e in edges:
			max_node_index += 1
			new_edge = (e[0], e[1], max_node_index)
			new_edges.append(new_edge)
			
			# Now actually add the nodes if necessary

			# Stokes 2nd order nodes should have no drag
			if stokes != None:
				stokes.add_node(0.0)
			
			if node != None:
				nnpos = 0.5 * (node.pos[e[0]] + node.pos[e[1]])
				node.add_node(nnpos)

		# Now, rebuild topology (and surface, maybe)
		# Make a dictionary
		new_edges = dict(((e[0],e[1]), e[2]) for e in new_edges)
		for i in range(self.num_elements):

			# Replace element with a second order one
			self.upgrade_element(i)

			# Get local edges and append nodes to element
			n = self.element[i].n
			#n = sorted(n)

			# 6 edges
			for j in range(3):
				for k in range(j + 1, 4):
					ed = (n[j], n[k])
					try:
						self.element[i].n.append(new_edges[ed])
					except(KeyError):
						ed = (n[k], n[j])
						self.element[i].n.append(new_edges[ed])
		# Now surface
		if surf != None:
			for i in range(surf.num_faces):

				# Replace face with a second order one, then convert to 4 linear ones
				surf.upgrade_face(i)

				# Get local edges and append nodes to element
				n = surf.face[i].n
				#n = sorted(n)

				# 3 edges
				for j in range(2):
					for k in range(j + 1, 3):
						ed = (n[j], n[k])
						try:
							surf.face[i].n.append(new_edges[ed])
						except(KeyError):
							ed = (n[k], n[j])
							surf.face[i].n.append(new_edges[ed])
			# Now make linear again
			for i in range(surf.num_faces):
				surf.split_face(0)			# Always zero, as split_face deletes the face at index

	def upgrade_element(self, index):

		# Replace element with a higher order one
		e = FFEA_element_tet_sec()
		e.n = self.element[index].n
		
		self.element.insert(index, e)
		self.element.pop(index + 1)

	#def check_normals(self, node, surf=None):
	#	
	#	# Get all faces on all elements and check that the normals are correctly oriented (pointing outwards)
	#	# Keep interior nodes interior!
#
#		num_elements_checked = 0
#		num_elements_correct = 0
#		while num_elements_correct != self.num_interior_elements:
#			num_elements_correct = 0
#			num_elements_checked = 0
#	
#			for i in range(self.num_interior_elements):
#				
#				# Get the 4 faces
#				efaces = [self.element[i].get_linear_face[n] for n in range(4)]
#
#				# For each face, check normal. If any one normal is wrong, they are all wrong. A switch of two indices will fix it.
#				# Trust me on this, or do some index permutations on a piece of paper :)

	def cull_interior(self, limitvol, node, surf=None):

		# First, get the interior stuff sorted, if we can
		self.calculateInterior(surf=surf)
		node.calculateInterior(top=self, surf=surf)

		# Now, for all interior elements, cull if vol < limitvol
		culled_elements = 0
		completed = False
		while(not completed):
			completed = True

			for i in range(self.num_elements):
				e = self.element[i]
				if not e.interior:
					continue
		
				vol = e.calc_volume(node)
				if vol < limitvol:
					
					#completed = False

					# Get the new node to add
					cent = e.calc_centroid(node)
					
					# Get the 4 adjacent elements we need to delete in addition to this one
					elstodelete = [i]
					for j in range(self.num_elements):
						#if e.sharesaface(self.element[j]):
						#	elstodelete.append(j)
						if e.sharesanedge(self.element[j]):
							elstodelete.append(j)

					nodestorenumber = sorted(e.n[0:4])

					# Delete elements (all interior)
					for j in reversed(elstodelete):
						self.element.pop(j)
					self.num_elements -= len(elstodelete)
					self.num_interior_elements -= len(elstodelete)
					culled_elements += len(elstodelete)

					# Delete all necessary nodes (not all interior :( )
					for j in reversed(nodestorenumber):
						if j < node.num_surface_nodes:
							node.num_surface_nodes -= 1
						else:
							node.num_interior_nodes -= 1
						
						try:
							node.pos.pop(j)
						except:
							node.pos = np.delete(node.pos, j, axis=0)

					node.num_nodes -= 4

					# Add new node (which will be interior!) and store new index
					node.add_node(cent)
					newnin = node.num_nodes - 1

					# Now renumber all nodes (to new node if in old element, else, to the new position in the node list)
					for j in range(self.num_elements):
						for k in range(len(self.element[j].n)):
							if self.element[j].n[k] in nodestorenumber:
								self.element[j].n[k] = newnin
							else:
								if self.element[j].n[k] < nodestorenumber[0]:
									pass
								elif self.element[j].n[k] < nodestorenumber[1]:
									self.element[j].n[k] -= 1
								elif self.element[j].n[k] < nodestorenumber[2]:
									self.element[j].n[k] -= 2
								elif self.element[j].n[k] < nodestorenumber[3]:
									self.element[j].n[k] -= 3
								else:
									self.element[j].n[k] -= 4

					break

					'''
					# Cull by replacing element with a node at the center
					# Add node (it will be interior)				
					cent = e.calc_centroid()
					node.add_node(cent, nodetype = 1)

					# Reconnect all old nodes to this node
					for j in range(self.num_elements):
						for k in range(4):
							if self.element[j].n[k] in e.n[0:4]:
								self.element[j].n[k] = node.num_nodes	# New index

					# Delete all old nodes and element (in reverse order so as to not ruin everything)
					for n in reversed(sorted(e.n)):
						node.pop(n)
						node.num_nodes -= 1
						node.num_interior_nodes -= 1

					self.element.remove(e)
					self.num_elements -= 1
					self.num_interior_elements -= 1
					break
					'''
		print ("Culled %d elements with volume < %e." % (culled_elements, limitvol))

	def get_smallest_lengthscale(self, node):

		length = float("inf")
		for e in self.element:
			new_length = e.get_smallest_lengthscale(node)
			if new_length < length:
				length = new_length

		return length

	def calculate_strain_energy(self, frame, frame0, mat):
		
		se = 0.0
		index = 0
		for el in self.element:
			se += el.calculate_strain_energy(frame, frame0, mat.element[index])
			index += 1
		return se

	def print_details(self):

		print ("num_elements = %d" % (self.num_elements))
		print ("num_surface_elements = %d" % (self.num_surface_elements))
		print ("num_interior_elements = %d" % (self.num_interior_elements))
		sleep(1)

		for e in self.element:
			index = self.element.index(e)
			outline = "Element " + str(index) + " "
			if(index < self.num_surface_elements):
				outline += "(Surface): "
			else:
				outline += "(Interior): "
			for n in e.n:
				outline += str(n) + " "

			print (outline)

	def write_to_file(self, fname):


		print("Writing to " + fname + "...")

		# Write differently depending on format
		base, ext = os.path.splitext(fname)

		if ext == ".vol":
			fout = open(fname, "a")
			fout.write("#  matnr      np      p1      p2      p3      p4\nvolumeelements\n%d\n" % (self.num_elements))
			for e in self.element:
				fout.write("1 %d" % (len(e.n)))
				for n in e.n:
					fout.write(" %d" % (n + 1))
				fout.write("\n")

			fout.write("\n\n")
		
		elif ext == ".top":
			fout = open(fname, "w")
			fout.write("ffea topology file\nnum_elements %d\nnum_surface_elements %d\nnum_interior_elements %d\n" % (self.num_elements, self.num_surface_elements, self.num_interior_elements))
			fout.write("surface elements:\n")
			for i in range(self.num_surface_elements):
				for n in self.element[i].n:
					fout.write("%d " % (n))
				fout.write("\n")
			fout.write("interior elements:\n")
			for i in range(self.num_surface_elements, self.num_elements):
				for n in self.element[i].n:
					fout.write("%d " % (n))
				fout.write("\n")
		else:
			print("Extension not recognised")
			raise IOError


		fout.close()
		print("done!")

	def calc_mass(self, mat, node, scale = 1.0):
	
		mass = 0.0
		index = 0
		for e in self.element:
			mass += e.calc_volume(node, scale) * mat.element[index][0]
			index += 1
		return mass

	# Takes index list of type intype ("node", "surf" etc) and returns the element list corresponding to those
	def index_switch(self, inindex, intype, limit=1, surf=None):
		
		outindex = []
		inindex = set(inindex)

		if intype.lower() == "node" or intype.lower() == "nodes":
			# Check if at least 'limit' nodes are in element
			for i in range(self.num_elements):
				if len(inindex & set(self.element[i].n)) >= limit:
		   			outindex.append(i)

		elif (intype.lower() == "surf" or intype.lower() == "surface" or intype.lower() == "face") and surf != None:
			
			# Check if face is on element
			outindex = [surf.face[i].elindex for i in inindex]

		else:
			raise IndexError

		return outindex

	def reset(self):

		self.CoM = None
		self.element = []
		self.num_elements = 0
		self.num_surface_elements = 0
		self.num_interior_elements = 0
		self.valid = False
		self.empty = True

class FFEA_element:

	def __init__(self):

		self.reset()

	def set_indices(self, alist):
		
		# Test for correct number of nodes
		if len(alist) != len(self.n):
			raise IndexError("Incorrect number of nodes for assignment to this element type.")

		for i in range(len(alist)):
			self.n[i] = int(alist[i])

	def calc_centroid(self, node):
	
		centroid = np.array([0.0,0.0,0.0])
		for i in self.n:
			centroid += node.pos[i]
			
		return centroid * (1.0 / len(self.n))
	
	def get_centroid(self):
		return self.centroid

	def get_linear_face(self, index, obj=True):
		
		# Define face i as the face that doesn't have node i in it
		if index == 0:
			n = [self.n[1], self.n[3], self.n[2]]
		elif index == 1:
			n = [self.n[0], self.n[2], self.n[3]]
		elif index == 2:
			n = [self.n[0], self.n[3], self.n[1]]
		elif index == 3:
			n = [self.n[0], self.n[1], self.n[2]]

		# Return either a face object, or a node list
		if obj:
			face = FFEA_surface.FFEA_face_tri_lin()
			face.set_indices(n)
			return face
		else:
			return n

	def sharesanedge(self, e):
	
		if len(set(e.n[0:4]) & set(self.n[0:4])) > 1:
			return True
		else:	
			return False

	def sharesaface(self, e):

		if len(set(e.n[0:4]) & set(self.n[0:4])) == 3:
			return True
		else:	
			return False

	def calc_volume(self, node, scale = 1.0):
		e = []
		for i in range(3):
			e.append(node.pos[self.n[i + 1]] - node.pos[self.n[0]])		

		return np.fabs(np.dot(e[2], np.cross(e[1], e[0])) / 6.0) * np.power(scale, 3.0)

	def calc_jacobian(self, node, scale = 1.0):
		J = np.array([[0.0 for i in range(3)] for j in range(3)])
		J[0][0] = node.pos[self.n[1]][0] - node.pos[self.n[0]][0]
		J[0][1] = node.pos[self.n[1]][1] - node.pos[self.n[0]][1]
		J[0][2] = node.pos[self.n[1]][2] - node.pos[self.n[0]][2]

		J[1][0] = node.pos[self.n[2]][0] - node.pos[self.n[0]][0]
		J[1][1] = node.pos[self.n[2]][1] - node.pos[self.n[0]][1]
		J[1][2] = node.pos[self.n[2]][2] - node.pos[self.n[0]][2]

		J[2][0] = node.pos[self.n[3]][0] - node.pos[self.n[0]][0]
		J[2][1] = node.pos[self.n[3]][1] - node.pos[self.n[0]][1]
		J[2][2] = node.pos[self.n[3]][2] - node.pos[self.n[0]][2]
		
		return J

	def calculate_strain_energy(self, frame, frame0, matel):

		# Initisalise (stuff we will need)
		se = 0.0
		startJ = self.calc_jacobian(frame0)
		J = self.calc_jacobian(frame)
		F = np.dot(J, np.linalg.inv(startJ))
		F = F.transpose()
		dF = np.linalg.det(F)
		G = matel[3]
		K = matel[4]
		C = K - (2.0/3.0) * G

		# Energy terms
		se += 0.5 * matel[3] * (np.dot(F,F.transpose()).trace() - 3)
		se += (C / 4.0) * (dF**2 - 1)
		se -= (0.5 * C + G) * np.log(dF)
		
		# Scale by volume
		se *= self.calc_volume(frame0)

		return se

	def get_smallest_lengthscale(self, node):

		# Smallest length is smalles node to opposite plane normal distance
		length = float("inf")
		for i in range(4):
			
			# Get a face
			f = self.get_linear_face(i)
			p = node.pos[self.n[i]]
			otherp = node.pos[self.n[(i + 1) % 4]] # A point in the plane (anything other than i in this element)
			n = f.calc_normal(node)

			# Define plane a plane[i]x_i + plane[3]
			plane = [n[0], n[1], n[2], -1 * np.dot(n, otherp)]
			
			# Distance is on the internet somewhere
			distance = np.fabs(np.dot(n, p) + plane[3])
			if distance < length:
				length = distance

		return length


	def reset(self):
		
		self.n = []
		self.interior = None

class FFEA_element_tet_lin(FFEA_element):

	def reset(self):

		self.n = [0,1,2,3]
		self.interior = None

class FFEA_element_tet_sec(FFEA_element):

	def reset(self):

		self.n = [0,1,2,3,4,5,6,7,8,9]
		self.interior = None
