import sys, os
from time import sleep
import numpy as np
from FFEA_exceptions import *

class FFEA_node:

	def __init__(self, fname = "", frame = 0):
	
		self.reset()

		# Empty fname give an empty object
		if fname == "":
			return

		try:
			self.load(fname)

		except FFEAFormatError as e:
			self.reset()
			print_error()
			print "Formatting error at line " + e.lin + "\nLine(s) should be formatted as follows:\n\n" + e.lstr
			raise

		except FFEAIOError as e:
			self.reset()
			print_error()
			print "Input error for file " + e.fname
			if e.fext != [""]:
				print "       Acceptable file types:"
				for ext in e.fext:
					print "       " + ext
		except IOError:
			raise

	def load(self, fname, findex = 0):

		sys.stdout.write("Loading FFEA node file...")
	
		# File format?
		base, ext = os.path.splitext(fname)
		try:
			if ext == ".node":

				# Check if tetgen
				with open(fname, "r") as fin:
					line = fin.readline().strip()
				if line == "ffea node file" or line == "walrus node file":
					self.load_FFEA_node(fname)
				elif len(line.split()) == 4:
					self.load_tetgen_node(fname)
				else:
					raise FFEAFormatError(lin=1)

			elif ext == ".out" or ext == ".traj" or ext == ".ftj":
				self.load_traj(fname, findex)

			elif ext == ".vol":
				self.load_vol(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".node", ".out", ".traj", ".ftj", ".vol"])

		except:
			raise

		self.valid = True
		sys.stdout.write("done!\n")

	def load_FFEA_node(self, fname):

		# Open file
		try:
			fin = open(fname, "r")

		except(IOError):
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea node file" and line != "walrus node file":
			raise FFEAFormatError(lin=1, lstr="ffea node file")

		try:
			num_nodes = int(fin.readline().split()[1])
			num_surface_nodes = int(fin.readline().split()[1])
			num_interior_nodes = int(fin.readline().split()[1])

		except IndexError:
			raise FFEAFormatError(lin="2-4", lstr="num_nodes %d\nnum_surface_nodes %d\nnum_interior_nodes %d")

		if fin.readline().strip() != "surface nodes:":
			raise FFEAFormatError(lin="5", lstr="surface nodes:")

		# Read nodes now
		try:
			j = 0
			for i in range(num_surface_nodes):
				sline = fin.readline().split()
				n = [float(sline[0]), float(sline[1]), float(sline[2])]
				self.add_node(n, nodetype = 0)

			if fin.readline().strip() != "interior nodes:":
				if num_interior_nodes != 0:
					raise FFEAFormatError(lin=num_surface_nodes + 6, lstr="interior nodes:")

			i = num_surface_nodes
			for j in range(num_interior_nodes):
				sline = fin.readline().split()
				n = [float(sline[0]), float(sline[1]), float(sline[2])]
				self.add_node(n, nodetype = 1)

		except (IndexError, ValueError):
			raise FFEAFormatError(lin=i+j+6, lstr="%f %f %f")
		except:
			raise

		fin.close()

		# Numpy it up, for speed
		self.pos = np.array(self.pos)

	def load_tetgen_node(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Test format
		sline = fin.readline().split()
		if int(sline[0]) == 1:
			pos = fin.tell()
			sline = fin.readline().split()
			try:
				if int(sline[0]) == 2:

					# Missing first line
					raise FFEAFormatError(lin=1, lstr="<num_nodes> <num_dimensions> 0 0")

			except(IndexError):
				# Missing first line
				raise FFEAFormatError(lin=1, lstr="<num_nodes> <num_dimensions> 0 0")

			fin.seek(pos)
			sline = fin.readline().split()
		try:
			num_nodes = int(sline[0])
			num_surface_nodes = 0
			num_interior_nodes = num_nodes

		except(IndexError, ValueError):
			raise FFEAFormatError(lin=1, lstr="<num_nodes> <num_dimensions> 0 0")

		# Read nodes now	
		for i in range(num_nodes):
			try:
				sline = fin.readline().split()
				if sline[0].strip() == "#":
					raise FFEAFormatError(lin=i + 2, lstr=str(i + 1) + " %f %f %f")

				sline = sline[1:]

				# Get a node
				n = [float(sline[0]), float(sline[1]), float(sline[2])]
				
				self.add_node(n)

			except(IndexError, ValueError):
				raise FFEAFormatError(lin=i + 2, lstr=str(i + 1) + " %f %f %f")


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
		try:
			if isinstance(self.pos, list):
				self.pos.append(n)
			else:
				self.pos = np.append(self.pos, [n], axis=0)
		except IndexError, ValueError:
			raise

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
	
	def write_to_file(self, fname, surf=None):

		print "Writing to " + fname + "..."

		# Write differently depending on format
		base, ext = os.path.splitext(fname)

		if ext == ".vol":
			fout = open(fname, "a")
			fout.write("#          X             Y             Z\npoints\n%d\n" % (self.num_nodes))
			for p in self.pos:
				fout.write("%22.16f  %22.16f  %22.16f\n" % (p[0], p[1], p[2]))

			fout.write("\n\n")

		elif ext == ".node":
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

		elif ext == ".obj":
			if surf == None:
				print "Error. Cannot write to '.obj' format without an associated 'surf' object"
				raise IOError
			else:
				surf.write_to_file(fname, node=self)
		else:
			print "Extension not recognised"
			raise IOError

		fout.close()
		print "done!"

	def scale(self, factor):
		
		self.pos *= factor

	def calc_centroid(self):
		
		return (1.0 / self.num_nodes) * np.sum(self.pos, axis = 0)
	
	def calc_CoM(self, top, mat):

		CoM = np.array([0.0,0.0,0.0])

		eindex = -1
		for e in top.element:
			eindex += 1
			elmass = e.calc_volume(self) * mat[eindex][0]
			CoM += elmass * (np.mean([self.pos[n] for n in e.n[0:4]]))
		self.CoM = CoM * 1.0/num_elements
		return self.CoM

	def get_CoM(self):
		return self.CoM

	def get_centroid(self):

		return self.centroid
	
	def translate(self, trans):
		self.pos += np.array(trans)
	
	def set_pos(self, pos):
		self.translate(np.array(pos) - self.calc_centroid())

	def rotate(self, rot):
		
		rot = np.array(rot)
		print rot
		
		# Translate to origin
		origin_trans = np.array([0.0,0.0,0.0]) - self.calc_centroid()
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
		self.centroid = None
		self.CoM = None
		self.num_nodes = 0
		self.num_surface_nodes = 0
		self.num_interior_nodes = 0
