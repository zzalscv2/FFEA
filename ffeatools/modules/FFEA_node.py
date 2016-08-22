import os
from time import sleep
import numpy as np

class FFEA_node:

	def __init__(self, fname = "", frame = 0):
	
		self.reset()
		if fname == "":
			return

		try:
			self.load(fname)
		except:
			return

	def load(self, fname, findex = 0):

		print("Loading FFEA node file...")

		# Test file exists
		if not os.path.exists(fname):
			print("\tFile '" + fname + "' not found. Returning empty object...")
			return
	
		# File format?
		base, ext = os.path.splitext(fname)
		if ext == ".node":

			# Check if tetgen
			try:
				with open(fname, "r") as fin:
					line = fin.readline().strip()
					if line == "ffea node file" or line == "walrus node file":
						
						try:
							self.load_FFEA_node(fname)
						except:
							print("\tUnable to load FFEA_node from " + fname + ". Returning empty object...")
					else:
						try:
							self.load_tetgen_node(fname)
						except:
							print("\tUnable to load FFEA_node from " + fname + ". Returning empty object...")
			except:
				print("\tUnable to load FFEA_node from " + fname + ". Returning empty object...")

		elif ext == ".out" or ext == ".traj":
			try:
				self.load_traj(fname, findex)
			except:
				print("\tUnable to load FFEA_node from " + fname + ", frame " + str(findex) + ". Returning empty object...")

		elif ext == ".vol":
			try:
				self.load_vol(fname)
			except:
				print("\tUnable to load FFEA_node from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_FFEA_node(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea node file" and line != "walrus node file":
			print("\tExpected 'ffea node file' but found " + line)
			raise TypeError

		num_nodes = int(fin.readline().split()[1])
		num_surface_nodes = int(fin.readline().split()[1])
		num_interior_nodes = int(fin.readline().split()[1])

		fin.readline()

		# Read nodes now
		nodetype = 0	
		while(True):
			sline = fin.readline().split()

			# Get a node
			try:
				if sline[0].strip() == "interior":
					nodetype = 1
					continue
	
				n = [float(sline[0]), float(sline[1]), float(sline[2])]

			except(IndexError):
				break

			self.add_node(n, nodetype = nodetype)

		fin.close()

		# Numpy it up, for speed
		self.pos = np.array(self.pos)

	def load_tetgen_node(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		sline = fin.readline().split()
		if len(sline) != 4:
			print("\tExpected '<num_nodes> <num_dimensions> 0 0' but found " + line)
			raise TypeError

		num_nodes = int(sline[0])
		num_surface_nodes = 0
		num_interior_nodes = num_nodes

		# Read nodes now	
		while(True):
			sline = fin.readline().split()

			if sline[0].strip() == "#":
				break

			sline = sline[1:]
			# Get a node
			try:
				n = [float(sline[0]), float(sline[1]), float(sline[2])]

			except(IndexError):
				break

			self.add_node(n)

		fin.close()

		# Numpy it up, for speed
		self.pos = np.array(self.pos)

	def load_vol(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		lines = fin.readlines()
		fin.close()

		# Test format and get to write place
		i = 0
		line = lines[i].strip()
		while line != "points":
			i += 1
			try:
				line = lines[i].strip()
			except(IndexError):
				print("\tCouldn't find 'points' line. File '" + fname + "' not formatted correctly.")
				self.reset()
				return

			continue

		# Get num_nodes
		i += 1
		num_nodes = int(lines[i])
		for j in range(i + 1, i + 1 + num_nodes):
			try:
				sline = lines[j].split()
				self.add_node([float(sline[0]), float(sline[1]), float(sline[2])])	
			except:
				print("\tCouldn't find the specified %d nodes. Only found %d. File '" + fname + "' not formatted correctly." % (num_nodes, i))
				self.reset()
				return

		# Numpy stuff up
		self.pos = np.array(self.pos)

	def add_node(self, n, nodetype = -1):

		# Numpy or not?
		if isinstance(self.pos, list):
			self.pos.append(n)
		else:
			self.pos = np.append(self.pos, [n], axis=0)

		self.num_nodes += 1
		
		if nodetype == -1:
			self.num_surface_nodes += 1
		elif nodetype == 0:
			self.num_surface_nodes += 1
		else:
			self.num_interior_nodes += 1

	def calculateInterior(self, top=None, surf=None):

		# We must have a topology and an associated surface, otherwise interior makes no sense
		if top == None or surf == None:
			print "Error. Cannot proceed without both a topology and a surface."
			return
		
		# Don't continue if we're already done
		if self.num_nodes == self.num_surface_nodes:
			return

		# Use surface to determine which nodes are interior and build a map
		amap = [-1 for i in range(self.num_nodes)]
		index = 0

		# Surface
		for f in surf.face:
			for n in f.n:
				if amap[n] == -1:
					amap[n] = index
					index += 1

		self.num_surface_nodes = index
		self.num_interior_nodes = self.num_nodes - self.num_surface_nodes

		# Now remainder are interior
		for i in range(self.num_nodes):
			if amap[i] == -1:
				amap[i] = index
				index += 1

		# Alter order of nodes
		oldpos = self.pos
 		self.pos = [np.array([0.0,0.0,0.0]) for i in range(self.num_nodes)]

		for n in range(len(amap)):
			self.pos[amap[n]] = oldpos[n]

		# And reassign surface and topologies
		for i in range(top.num_elements):

			# Map node indices
			for j in range(len(top.element[i].n)):
				top.element[i].n[j] = amap[top.element[i].n[j]]

		for i in range(surf.num_faces):
			for j in range(len(surf.face[i].n)):
				surf.face[i].n[j] = amap[surf.face[i].n[j]]
	
		# And make sure the interior node of surface elements if at the end of the list of the linear indices
		#for i in range(top.num_interior_elements, top.num_elements):
		#	for j in range(4):
		#		if top.element[i].n[j] < node.num_interior_nodes:
		#			break
		#	
		#	# j is the index of the interior node
		#	if j < 3:
		#		# Permute the first 3 indices
		#		while j < 4
		#	else:
		#		# Permute the last 3 indices 

	def calculate_dimensions(self):
		
		# min, max
		dims = [[float("inf"), -1 * float("inf")] for i in range(3)]

		for p in self.pos:
			for i in range(3):
				if p[i] < dims[i][0]:
					dims[i][0] = p[i]
				if p[i] > dims[i][1]:
					dims[i][1] = p[i]

		return [d[1] - d[0] for d in dims]

	def print_details(self):

		print "num_nodes = %d" % (self.num_nodes)
		print "num_surface_nodes = %d" % (self.num_surface_nodes)
		print "num_interior_nodes = %d" % (self.num_interior_nodes)
		sleep(1)

		index = -1
		for n in self.pos:
			index += 1
			outline = "Node " + str(index) + " "
			if(index < self.num_surface_nodes):
				outline += "(Surface): "
			else:
				outline += "(Interior): "
			for xyz in n:
				outline += "%6.3e " % (xyz)
			
			print outline
	
	def write_to_file(self, fname):

		print "Writing to " + fname + "..."

		# Write differently depending on format
		base, ext = os.path.splitext(fname)

		if ext == ".vol":
			fout = open(fname, "a")
			fout.write("#          X             Y             Z\npoints\n%d\n" % (self.num_nodes))
			for p in self.pos:
				fout.write("%22.16f  %22.16f  %22.16f\n" % (p[0], p[1], p[2]))

			fout.write("\n\n")
		else:
			fout = open(fname, "w")
			fout.write("ffea node file\nnum_nodes %d\nnum_surface_nodes %d\nnum_interior_nodes %d\n" % (self.num_nodes, self.num_surface_nodes, self.num_interior_nodes))
		
			# Surface nodes
			fout.write("surface nodes:\n")
			for i in range(self.num_surface_nodes):
				fout.write("%6.3f %6.3f %6.3f\n" % (self.pos[i][0], self.pos[i][1], self.pos[i][2]))

			# Interior nodes
			fout.write("interior nodes:\n")
			for i in range(self.num_surface_nodes, self.num_nodes, 1):
				fout.write("%6.3f %6.3f %6.3f\n" % (self.pos[i][0], self.pos[i][1], self.pos[i][2]))
		fout.close()
		print "done!"

	def scale(self, factor):
		
		self.pos *= factor

	def get_centroid(self):
		
		return (1.0 / self.num_nodes) * np.sum(self.pos, axis = 0)
		
	def translate(self, trans):
		self.pos += np.array(trans)
	
	def set_pos(self, pos):
		self.translate(np.array(pos) - self.get_centroid())

	def rotate(self, rot):
		
		rot = np.array(rot)
		print rot
		
		# Translate to origin
		origin_trans = np.array([0.0,0.0,0.0]) - self.get_centroid()
		self.translate(origin_trans)
		
		if rot.size == 3:
		
			# Rotate in x, then y, then z
			c = np.cos
			s = np.sin
			x = np.radians(rot[0])
			y = np.radians(rot[1])
			z = np.radians(rot[2])
			Rx = np.array([[1, 0, 0],[0,c(x),-s(x)],[0,s(x),c(x)]])
			Ry = np.array([[c(y), 0, s(y)],[0,1,0],[-s(y),0,c(y)]])
			Rz = np.array([[c(z),-s(z),0],[s(z),c(z),0], [0,0,1]])
		
			# x, y, z. Change if you want
			R = np.dot(Rz, np.dot(Ry, Rx))
			
		elif rot.size == 9:
			R = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
			for i in range(3):
				for j in range(3):
					R[i][j] = rot[3 * i + j]
		
		else:
			return
						
		for i in range(self.num_nodes):
			self.pos[i] = np.dot(R, self.pos[i])
			
		# Translate back
		self.translate(-1 * origin_trans)

	def reset(self):

		self.pos = []
		self.num_nodes = 0
		self.num_surface_nodes = 0
		self.num_interior_nodes = 0
