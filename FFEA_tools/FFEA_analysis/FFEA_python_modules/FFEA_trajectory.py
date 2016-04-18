from os import path
import numpy as np
import FFEA_frame, FFEA_pdb
import sys

class FFEA_trajectory:

	def __init__(self, fname="", surf=None, load_all=1, frame_rate = 1):

		self.reset()

		# Return empty object if fname not initialised
		if fname == "":
			return

		if self.load(fname, load_all=load_all, surf=surf, frame_rate = frame_rate) == 1:
			print("\tLoading of '" + fname + "' failed. Returning empty object...")
			return	
		
	def load(self, fname, surf=None, load_all=1, frame_rate = 1):

		print("Loading FFEA trajectory file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found. Returning empty object...")
			return 1
	
		# Clear everything for beginning
		self.reset()

		# Header first, for sure
		if self.load_header(fname) == 1:
			print("\tUnable to load header information from FFEA_trajectory '" + fname + "'.")
			self.reset()
			return 1

		# Then rest of trajectory.
		if(load_all == 1):
			all_frames = -1
			while(True):
				all_frames += 1
				if all_frames % frame_rate != 0:
					if self.skip_frame() == -1:
						break
				
				elif(self.load_frame(surf=surf) != 0):
					print("done! Successfully read " + str(self.num_frames) + " frame/s from '" + fname + "'.")
					break

				if all_frames % 100 == 0:
					print "Frames parsed = ", str(all_frames)

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

				sline = self.traj.readline().split()
				cindex = int(sline[3][0])
				step = int(sline[5])

			except(IndexError):
				return 1

			except(ValueError):

				# Don't need to reset though!
				print("Unable to read conformation index for blob " + str(bindex) + " at frame " + str(self.num_frames))
				return 1
				
			# Get a frame first (STATIC should just continue, single frame can be added prior or after)
			if self.traj.readline().strip() == "DYNAMIC":
				b[cindex].motion_state = "DYNAMIC"
				frame = FFEA_frame.FFEA_frame()
			else:
				b[cindex].motion_state = "STATIC"
				continue

			# Read stuff
			frame.load_from_traj(self.traj)
			frame.set_step(step)

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

	def skip_frame(self):

		num_asterisks = 0
		while(True):
			line = self.traj.readline()
			try:
				if line[0] == "*":
					num_asterisks += 1
					if num_asterisks == 2:
						break
			except:
				return -1

	def build_from_pdb(self, pdb, scale = 1):

		# Single blob single conf
		self.set_header(1, [1], [[pdb.blob[0].num_atoms]])
		print scale
		for i in range(pdb.num_frames):
			frame = FFEA_frame.FFEA_frame()
			frame.pos = pdb.blob[0].frame[i].pos * scale
			frame.set_step(i)
			self.blob[0][0].num_nodes = len(frame.pos)
			self.blob[0][0].frame.append(frame)
			self.num_frames += 1

	# Manually set a single frame
	def set_single_frame(self, node, surf = None, step = 0):

		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
	
				# Reset all first
				self.blob[i][j].reset()

				# Add frame only for conf 0
				if j == 0:
					frame = FFEA_frame.FFEA_frame()
					frame.set_step(step)
					frame.pos = node[i].pos

					if surf != None:
						frame.calc_normals(surf[i])

					self.blob[i][j].num_nodes = len(frame.pos)
					self.blob[i][j].frame.append(frame)

				else:
					self.blob[i][j].frame.append(None)
		
		self.num_frames = 1
	
	def write_to_file(self, fname, frames=None, frame_rate = 1):

		# Get a file object
		fout = open(fname, "w")
		
		# Write header info
		self.write_header_to_file(fout)

		# Write frames
		if frames == None:
			frames = [0,self.num_frames]

		for i in range(frames[0], frames[1], frame_rate):
			self.write_frame_to_file(fout, i)

		fout.close()

	def write_header_to_file(self, fout):

		
		fout.write("FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs %d\nNumber of Conformations" % (self.num_blobs))
		for i in self.num_conformations:
			fout.write(" %d" % (i))
		fout.write("\n")
		for i in range(self.num_blobs):
			fout.write("Blob %d:" % (i))
			for j in range(self.num_conformations[i]):
				fout.write(" Conformation %d Nodes %d" % (j, self.num_nodes[i][j]))
			fout.write("\n")
		fout.write("\n*\n")

	def write_frame_to_file(self, fout, index):

		if index >= self.num_frames:
			return

		# Traj data
		cur_conf = []
		for b in self.blob:
			bindex = self.blob.index(b)
			for c in b:
				cindex = b.index(c)
				if c.motion_state == "STATIC":
					fout.write("Blob %d, Conformation %d, step %d\nSTATIC\n" % (bindex, cindex, 0))
					cur_conf.append(0)
					break
				else:
					if c.frame[index] == None:
						continue
				
					cur_conf.append(cindex)
					fout.write("Blob %d, Conformation %d, step %d\n" % (bindex, cindex, c.frame[index].step))
					fout.write(c.motion_state + "\n")
					c.frame[index].write_to_traj(fout)

		# Kinetic Data
		fout.write("*\nConformation Changes:\n")
		for b in self.blob:
			bindex = self.blob.index(b)
			try:
				for c in b:
					cindex = b.index(c)
					if c.frame[index + 1] != None:
						next_conf = cindex
						break
			except:
				next_conf = cur_conf[bindex]

			fout.write("Blob %d: Conformation %d -> Conformation %d\n" % (bindex, cur_conf[bindex], next_conf))
		fout.write("*\n")		

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
	
	def set_subblob(self, pin):

		if max(pin.index) >= self.num_nodes:
			print("Error. Pinned node index %d is larger than num_nodes, %d." % (max(pin.index), self.num_nodes))
			return None 
			
		self.subblob.append(pin.index)
		self.num_subblobs += 1
	
	def get_centroid_trajectory(self, subblob_index = -1):
		
		if subblob_index == -1:
			indices = [i for i in range(self.num_nodes)]
		else:
			try:
				indices = self.subblob[subblob_index]
			except(IndexError):
				print("Error. Subblob index %d out of range (num_subblobs = %d)." % (subblob_index, self.num_subblobs))
				return None, None

		# Build the trajectory
		subblob_size = len(indices)		
		ctraj = []
		step = []
		
		print len(self.frame)
		for f in self.frame:
			centroid = np.array([0.0,0.0,0.0])
			for i in indices:
				centroid += f.pos[i]

			centroid /= subblob_size
			ctraj.append(centroid)
			step.append(f.step)
			
		return np.array(step), np.array(ctraj)

	def reset(self):

		self.motion_state = "DYNAMIC"
		self.num_nodes = 0
		self.num_subblobs = 0
		self.frame = []
		self.subblob = []
