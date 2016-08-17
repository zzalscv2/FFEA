import os, sys
from time import sleep
import numpy as np
import FFEA_surface

class FFEA_topology:

	def __init__(self, fname = ""):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname):

		print("Loading FFEA topology file...")

		# Test file exists
		if not os.path.exists(fname):
			print("\tFile '" + fname + "' not found. Returning empty object...")
			return
	
		# File format?
		base, ext = os.path.splitext(fname)
		if ext == ".top":
			try:
				self.load_top(fname)
			except:
				print("\tUnable to load FFEA_topology from " + fname + ". Returning empty object...")

		elif ext == ".ele":
			try:
				self.load_ele(fname)
			except:
				print("\tUnable to load FFEA_topology from " + fname + ". Returning empty object...")

		elif ext == ".vol":
			try:
				self.load_vol(fname)
			except:
				print("\tUnable to load FFEA_topology from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_top(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
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
			self.add_element(el)

		fin.close()

	def load_vol(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
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
		for i in range(num_elements):
			sline = fin.readline().split()[2:]
			el = FFEA_element_tet_lin()
			el.set_indices(sline)
			self.add_element(el)

		# Indexing from 0
		for i in range(self.num_elements):
			for j in range(4):
				self.element[i].n[j] -= 1

	def load_ele(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
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
	
	def extract_surface(self):
		
		faces = []
		facessorted = []
		elindex = []
		surf = FFEA_surface.FFEA_surface()

		# Get all faces in entire system. If only occurs once, it's a surface face. Keep the order
		for i in range(self.num_elements):
			e = self.element[i]

			if len(e.n) == 10:
				print "Error. Functionality not yet available for 2nd order elements"
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
			print "Cannot calculate interior elements without an associated surface (currently)."
			return

		# Don't continue if we're already done
		for e in self.element:
			if e != None:
				return

		# Set all elements as default to interior elements
		self.num_interior_elements = self.num_elements
		self.num_surface_elements = 0

		for i in range(self.num_elements):
			self.element[i].interior = True

		# Get a map, so we know what element wil go where
		amap = [-1 for i in range(self.num_elements)]
		index = 0

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
			print "Element ", index, " does not exist."
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
				print "Error. Increasing to order > 2 is currently not supported"
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
		calculateInterior(self, surf=surf)

		# Now, for all interior elements, cull if vol < limitvol
		while(True):
			completed = True

			for i in range(self.num_elements):
				e = self.element
				if not e.interior:
					continue
		
				vol = e.calc_volume(node)
			
				if vol < limitvol:
					
					completed = False

					# Get the new node to add and it's index
					cent = e.calc_centroid()
					newnin = node.num_nodes

					# Get the 4 adjacent elements we need to delete
					elstodelete = []
					for j in range(self.num_elements):
						if e.sharesface(self.element[j]):
							elstodelete.append(j)

					
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
			if completed:
				break

	def print_details(self):

		print "num_elements = %d" % (self.num_elements)
		print "num_surface_elements = %d" % (self.num_surface_elements)
		print "num_interior_elements = %d" % (self.num_interior_elements)
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

			print outline

	def write_to_file(self, fname):


		print "Writing to " + fname + "..."

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
			print "Could not write topology to " + fname
			return

		fout.close()
		print "done!"

	def calc_mass(self, mat, node, scale = 1.0):
	
		mass = 0.0
		index = 0
		for e in self.element:
			mass += e.get_volume(node, scale) * mat.element[index][0]
			index += 1
		return mass

	def reset(self):

		self.element = []
		self.num_elements = 0
		self.num_surface_elements = 0
		self.num_interior_elements = 0

class FFEA_element:

	def __init__(self):

		self.reset()

	def set_indices(self, alist):
		
		# Test for correct number of nodes
		if len(alist) != len(self.n):
			print "Incorrect number of nodes for assignment to this element type."
			return

		for i in range(len(alist)):
			self.n[i] = int(alist[i])

	def calc_centroid(self, node):
	
		centroid = np.array([0.0,0.0,0.0])
		for i in self.n:
			centroid += node.pos[i]
			
		return centroid * (1.0 / len(self.n))
	
	def get_linear_face(self, index, obj=True):
		
		if index == 0:
			n = [self.n[0], self.n[1], self.n[2]]
		elif index == 1:
			n = [self.n[0], self.n[3], self.n[1]]
		elif index == 2:
			n = [self.n[0], self.n[2], self.n[3]]
		elif index == 3:
			n = [self.n[1], self.n[3], self.n[2]]

		# Return either a face object, or a node list
		if obj:
			face = FFEA_surface.FFEA_face_tri_lin()
			face.set_indices(n)
			return face
		else:
			return n

	def sharesface(self, index):
		return
	def calc_volume(self, node, scale = 1.0):
		e = []
		for i in range(3):
			e.append(node.pos[self.n[i + 1]] - node.pos[self.n[0]])		

		return np.fabs(np.dot(e[2], np.cross(e[1], e[0])) / 6.0) * np.power(scale, 3.0)

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
