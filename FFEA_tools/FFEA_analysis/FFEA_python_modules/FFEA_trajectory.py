from os import path
import numpy as np
import FFEA_frame

class FFEA_trajectory:

	def __init__(self, fname="", surf=None, load_all=1):

		self.reset()

		# Return empty object if fname not initialised
		if fname == "":
			return

		if self.load(fname, load_all=load_all, surf=surf) == 1:
			print("\tLoading of '" + fname + "' failed. Returning empty object...")
			return	
		
	def load(self, fname, surf=None, load_all=1):

		print("Loading FFEA trajectory file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found. Returning empty object...")
			return
	
		# Clear everything for beginning
		self.reset()

		# Header first, for sure
		if self.load_header(fname) == 1:
			print("\tUnable to load header information from FFEA_trajectory '" + fname + "'.")
			self.reset()
			return 1

		# Then rest of trajectory.
		if(load_all == 1):
			while(True):
				if(self.load_frame(surf=surf) != 0):
					print("done! Successfully read " + str(self.num_frames) + " frame/s from '" + fname + "'.")
					break

	def load_header(self, fname):

		# Get a file object and store it
		try:
			self.traj = open(fname, "r")

		except(IOError):
			print("\tFailed to open '" + fname + "' for reading.")
			return 1

		# Now, read only the information from the top of the file

		# Title
		line = self.traj.readline().strip()
		if line != "FFEA_trajectory_file":
			print("\tExpected to read 'FFEA_trajectory_file' but read '" + line + "'. This may not be an FFEA trajectory file.")
			self.reset()
			return 1


		self.traj.readline()
		self.traj.readline()

		# num_blobs
		try:
			line = self.traj.readline()
			self.num_blobs = int(line.split()[3])

		except(IndexError, ValueError):
			print("\tExpected to read 'Number of Blobs %d' but read '" + line + "'.")
			self.reset()
			return 1

		# num_conformations
		try:
			line = self.traj.readline()
			sline = line.split()[3:]
			self.num_conformations = [int(s) for s in sline]

		except(IndexError, ValueError):
			print("\tExpected to read 'Number of Conformations %d %d ....%d' but read '" + line + "'.")
			self.reset()
			return 1

		# num_nodes
		self.num_nodes = [[0 for i in range(self.num_conformations[j])] for j in range(self.num_blobs)]
		for i in range(self.num_blobs):
			try:
				line = self.traj.readline()
				sline = line.split()[2:]

				for j in range(self.num_conformations[i]):
					self.num_nodes[i][j] = int(sline[4 * j + 3])

			except(IndexError, ValueError):
				print("\tExpected to read 'Blob " + str(i) + ": Conformation 0 Nodes %d Conformation 1 Nodes %d....Conformation " + str(num_conformations[i] - 1) + " Nodes %d' but read '" + line + "'.")
				self.reset()
				return 1

		# final whitespace until '*'
		while(self.traj.readline().strip() != "*"):
			pass

		# Finally, build the objects
		self.blob = [[FFEA_traj_blob(self.num_nodes[i][j]) for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

	# Manually set header data
	def set_header(self, num_blobs, num_conformations, num_nodes):

		self.num_blobs = num_blobs
		self.num_conformations = num_conformations
		self.num_nodes = num_nodes

		# Still build the objects
		self.blob = [[FFEA_traj_blob(self.num_nodes[i][j]) for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

	# This function must be run as fast as possible! Error checking will be at a minimum. This function is standalone so it can be threaded
	def load_frame(self, surf=None):
	
		# For each blob
		for b in self.blob:
			try:
				# Get indices
				bindex = self.blob.index(b)
				cindex = int(self.traj.readline().split()[3][0])

			except(IndexError):
				return 1

			except(ValueError):

				# Don't need to reset though!
				print("Unable to read conformation index for blob " + str(bindex) + " at frame " + str(self.num_frames))
				return 1
				
			# Get a frame first (STATIC should just continue, single frame can be added prior or after)
			if self.traj.readline().strip() == "DYNAMIC":
				frame = FFEA_frame.FFEA_frame()
			else:
				continue

			# Read stuff
			frame.load_from_traj(self.traj)

			# Load normals if necessary
			#if surf != None:
			#	frame.calc_normals(surf[bindex][cindex])

			# Append frame
			b[cindex].frame.append(frame)

			# Append None to all frames that aren't active
			for c in range(self.num_conformations[bindex]):
				if c != cindex:
					b[c].frame.append(None)

		# Gloss over kinetics stuff
		self.traj.readline()
		while(self.traj.readline().strip() != "*"):
			pass

		self.num_frames += 1
		return 0

	# Manually set a single frame
	def set_single_frame(self, node, surf = None):

		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
	
				# Reset all first
				self.blob[i][j].reset()

				# Add frame only for conf 0
				if j == 0:
					frame = FFEA_frame.FFEA_frame()
					frame.pos = node[i].pos

					if surf != None:
						frame.calc_normals(surf[i])

					self.blob[i][j].num_nodes = len(frame.pos)
					self.blob[i][j].frame.append(frame)

				else:
					self.blob[i][j].frame.append(None)
		
		self.num_frames = 1
		
	def reset(self):

		self.traj = None
		self.num_frames = 0
		self.num_blobs = 0
		self.num_conformations = []
		self.num_nodes = []
		self.blob = []

class FFEA_traj_blob:

	def __init__(self, num_nodes=0):

		self.reset()
		self.num_nodes = num_nodes

	# Manually set a frame
	def set_frame(self, frame):
		self.frame.append(frame)
		
	def reset(self):

		self.num_nodes = 0
		self.frame = []
