from os import path
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
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".node":
			try:
				self.load_node(fname)
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

	def load_node(self, fname):

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

	def add_node(self, n, nodetype = 0):

		self.pos.append(n)
		self.num_nodes += 1
		
		if nodetype == 0:
			self.num_surface_nodes += 1
		else:
			self.num_interior_nodes += 1
		
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
