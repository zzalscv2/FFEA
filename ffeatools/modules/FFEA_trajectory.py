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

		self.load(fname, load_all=load_all, surf=surf, frame_rate = frame_rate, num_frames_to_read = num_frames_to_read, start = start)

		return	
		
	def load(self, fname, surf=None, load_all=1, frame_rate = 1, num_frames_to_read = 1000000, start = 0):

		print("Loading FFEA trajectory file...")

		# Test file exists
		if not path.exists(fname):
			raise IOError("No trajectory found at that location")
	
		# Clear everything for beginning
		self.reset()

		# Header first, for sure
		self.load_header(fname)

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
			raise IOError("\tFailed to open '" + fname + "' for reading.")

		# Now, read only the information from the top of the file

		# Title
		line = self.traj.readline().strip()
		if line != "FFEA_trajectory_file":
			raise IOError("\tExpected to read 'FFEA_trajectory_file' but read '" + line + "'. This may not be an FFEA trajectory file.")

		self.traj.readline()
		self.traj.readline()

		# num_blobs
		try:
			line = self.traj.readline()
			self.num_blobs = int(line.split()[3])

		except(IndexError, ValueError):
			raise IOError("\tExpected to read 'Number of Blobs %d' but read '" + line + "'.")

		# num_conformations
		try:
			line = self.traj.readline()
			sline = line.split()[3:]
			self.num_conformations = [int(s) for s in sline]

		except(IndexError, ValueError):
			raise IOError("\tExpected to read 'Number of Conformations %d %d ....%d' but read '" + line + "'.")

		# num_nodes
		self.num_nodes = [[0 for i in range(self.num_conformations[j])] for j in range(self.num_blobs)]
		for i in range(self.num_blobs):
			try:
				line = self.traj.readline()
				sline = line.split()[2:]

				for j in range(self.num_conformations[i]):
					self.num_nodes[i][j] = int(sline[4 * j + 3])

			except(IndexError, ValueError):
				raise IOError("\tExpected to read 'Blob " + str(i) + ": Conformation 0 Nodes %d Conformation 1 Nodes %d....Conformation " + str(num_conformations[i] - 1) + " Nodes %d' but read '" + line + "'.")

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
    
    	# ye who enter here: do not 'fix' this! The script is not handling an
	# exception poorly, it is asking for forgiveness, not permission.
	# Because it's faster.
	# Signed, someone who tried to 'fix' this.

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
		try:
			for i in range(self.num_blobs):
				for j in range(self.num_conformations[i]):
					del self.blob[i][j].frame[index]
				
			self.num_frames -= 1
		except:
			raise
			
	def build_from_pdb(self, pdb, scale = 1):

		# Single blob single conf
		self.set_header(1, [1], [[pdb.chain[0].num_atoms]])

		# Get all frames in one go (they're the same core objects) 
		self.blob[0][0].num_nodes = pdb.chain[0].num_atoms
		self.blob[0][0].frame = pdb.chain[0].frame
		
		index = 0
		for f in self.blob[0][0].frame:
			f.pos *= scale
			f.set_step(index)
			index += 1

		self.num_frames = pdb.chain[0].num_frames

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

		print("Writing trajectory to file\n\tData will be written to %s\n" % (fname))

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
		bindex = 0
		for b in self.blob:
			next_conf = cur_conf[bindex]
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
			bindex += 1
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
		"""
		Create a subblob from a pin object.
		In: self, the pin object
		Out: the ID of the generated subblob.
		"""

		if max(pin.index) >= self.num_nodes:
			raise IndexError("Error. Pinned node index %d is larger than num_nodes, %d." % (max(pin.index), self.num_nodes))

		self.subblob.append(pin.index)
		self.num_subblobs += 1
		return self.num_subblobs # let the user assign a nice friendly name to their subblob
	
	def calc_centroid_trajectory(self, subblob_index = -1):
		"""
		Calculate the centroid (average position of all the nodes) of a
		given blob. You can also specify the index of a sub-blob, and get
		the centroid of that sub-blob.
		In: self, subblob index
		Out: a 2-d array, which is [x,y,z] wide and as tall as the number
		of nodes in your blob\sub-blob.
		"""
		
		if subblob_index == -1:
			indices = [i for i in range(self.num_nodes)]
		else:
			try:
				indices = self.subblob[subblob_index]
			except(IndexError):
				raise IndexError("Error. Subblob index %d out of range (num_subblobs = %d)." % (subblob_index, self.num_subblobs))

		# Build the trajectory
		subblob_size = len(indices)		
		ctraj = []
		step = []
		
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

