from os import path
import numpy as np
import FFEA_frame, FFEA_pdb
import sys

class FFEA_trajectory:

	def __init__(self, fname="", surf=None, load_all=1, frame_rate = 1, num_frames_to_read = 1000000, start = 0):

		self.reset()

		# Return empty object if fname not initialised
		if fname == "" or fname == None:
			return

		if self.load(fname, load_all=load_all, surf=surf, frame_rate = frame_rate, num_frames_to_read = num_frames_to_read, start = start) == 1:
			print("\tLoading of '" + fname + "' failed. Returning empty object...")
		else:
			self.valid = True

		return	
		
	def load(self, fname, surf=None, load_all=1, frame_rate = 1, num_frames_to_read = 1000000, start = 0):

		print("Loading FFEA trajectory file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found. Returning empty object...")
			self.reset()
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
			all_frames = 0
			while(True):

				# Have we read enough frames?
				#print all_frames, num_frames_to_read
				if((all_frames - start) == num_frames_to_read):
					print("\ndone! Successfully read " + str(self.num_frames) + " frame/s from '" + fname + "'.")
					break

				# Skip or load
				if (all_frames - start) % frame_rate != 0 or all_frames < start:
					if self.skip_frame() == 1:
						print("\ndone! Successfully read " + str(self.num_frames) + " frame/s from '" + fname + "'.")
						break
				
				elif(self.load_frame(surf=surf) != 0):
					print("\ndone! Successfully read " + str(self.num_frames) + " frame/s from '" + fname + "'.")
					break

				all_frames += 1
				#if self.num_frames % 100 == 0:
				#	print "Frames parsed = ", str(all_frames)

				sys.stdout.write("\rFrames read = %d, Frames skipped = %d" % (self.num_frames, all_frames - self.num_frames))
				sys.stdout.flush()

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

		# final whitespace until '*' and save the file pos
		while(self.traj.readline().strip() != "*"):
			pass

		self.fpos = self.traj.tell()
		
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
		eof = False
		bindex = -1
		for b in self.blob:
			try:
				# Get indices
				bindex += 1

				sline = self.traj.readline().split()
				cindex = int(sline[3][0])
				step = int(sline[5])

			except(IndexError):
				self.traj.seek(self.fpos)
				return 1

			except(ValueError):

				# Don't need to reset though!
				self.traj.seek(self.fpos)
				print("Unable to read conformation index for blob " + str(bindex) + " at frame " + str(self.num_frames))
				return 1
				
			# Get a motion_state
			b[cindex].motion_state = self.traj.readline().strip()
			#if self.traj.readline().strip() == "DYNAMIC":
		#		b[cindex].motion_state = "DYNAMIC"
		#		frame = FFEA_frame.FFEA_frame()
		#	else:
		#		b[cindex].motion_state = "STATIC"
		#		frame = b[cindex].frame[0]
		#		continue

			# Do different things depending on motion state recieved
			if b[cindex].motion_state == "STATIC":
				
				# Nothing to read; frame cannot be read from trajectory
				frame = None

			else:

				# Get a frame
				frame = FFEA_frame.FFEA_frame()

				# Try to read stuff
				success = frame.load_from_traj(self.traj)
			
				# We are at eof, or halfway through a frame being written
				if success == 1:
					eof = True
					break
				
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

		if not eof:
		
			# Gloss over kinetics stuff
			self.traj.readline()
			while(self.traj.readline().strip() != "*"):
				pass

			self.num_frames += 1
			self.fpos = self.traj.tell()
			return 0
		else:
			self.traj.seek(self.fpos)

	def scale(self, factor, frame_index):

		for b in range(self.num_blobs):
			for c in range(self.num_conformations[b]):
				try:
					self.blob[b][c].frame[frame_index].scale(factor)
				except:
					continue	

	def translate(self, trans):
		for b in range(self.num_blobs):
			for c in range(self.num_conformations[b]):
				for f in range(self.num_frames):
					self.blob[b][c].frame[f].translate(trans)

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
				return 1

		return 0

	def delete_frame(self, index=-1):
		
		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
				del self.blob[i][j].frame[index]
				
		self.num_frames -= 1
				
	def build_from_pdb(self, pdb, scale = 1):

		# Single blob single conf
		self.set_header(1, [1], [[pdb.blob[0].num_atoms]])
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

		print "Writing trajectory to file\n\tData will be written to %s\n" % (fname)

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
		bindex = -1
		for b in self.blob:
			bindex += 1
			cindex = -1
			for c in b:
				cindex += 1
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
		bindex = -1
		for b in self.blob:
			bindex += 1
			cindex = -1
			try:
				for c in b:
					cindex += 1
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
		self.valid = False

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


# External functions
def get_num_frames(fname):

	fin = open(fname, "r")
	if fin.readline().strip() != "FFEA_trajectory_file":
		print("\tExpected to read 'FFEA_trajectory_file' but read '" + line + "'. This may not be an FFEA trajectory file.")
		self.reset()
		return 1
		
	num_asterisks = 0	
	for line in fin:
		if "*" in line:
			num_asterisks += 1
	
	return (num_asterisks - 1) / 2

