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

import sys, os
from time import sleep
import numpy as np
from FFEA_exceptions import *
import FFEA_frame

class FFEA_pdb:

	def __init__(self, fname = "", num_frames = 100000):
	
		self.reset()

		if fname == "":
			self.valid = True
			sys.stdout.write("Empty pdb object initialised.\n")
			return

		try:
			self.load(fname, num_frames = num_frames)

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

	def load(self, fname, num_frames = 100000):

		sys.stdout.write("Loading FFEA pdb file...\n")

		# File format?
		base, ext = os.path.splitext(fname)
		try:
			if (ext == ".pdb"):
				self.load_pdb(fname, num_frames_to_read = num_frames)
			else:
				raise FFEAIOError(fname=fname, fext=[".pdb"])

		except:
			raise
		
		self.valid = True
		self.empty = False
		sys.stdout.write("...done!\n")

	def load_pdb(self, fname, num_frames_to_read = 100000):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise IOError

		#
		# Try to make robust against dodgy formating by the user
		#
		sys.stdout.write("\tScanning first frame to determine format...")
		
		# Read to first bit of interest (each MODEL is a frame, each chain a blob)
		line = fin.readline()
		while True:

			# Ignore MODEL, only insert when we output
			if line[0:4] == "ATOM":
				aID=int(line[6:12])
				cID = line[21]
				fin.seek(-len(line),1)
				start_pos = fin.tell()
				break
			else:
				line = fin.readline()
		num_atoms = [0]
		num_chains = 1

		# Now read the first frame (can't trust TERs to exist, sometimes people leave them out)
		while True:
			line = fin.readline()
			if line[0:4] == "ATOM":

				# New frame?
				if int(line[6:12]) < aID:
					break
				
				aID = int(line[6:12])

				# New chain?
				if line[21] != cID:
					cID = line[21]
					num_chains += 1
					num_atoms.append(0)

				num_atoms[-1] += 1

			elif line[0:3] == "TER":

				# See if next line gives a new chain. If not, register a new chain
				check_pos = fin.tell()
				line = fin.readline()
				
				## Any failure to read next line constitues the end of the frame
				try:
					if line[21] == cID:
						num_chains += 1
						num_atoms.append(0)
				except:
					break

				fin.seek(check_pos)

			elif line[0:3] == "END" or line.strip() == "":
				break

		sys.stdout.write("\tdone!\n")

		#
		# Reset to beginning and read atomic structures
		#
		sys.stdout.write("\tReading first frame to determine structure...")		
		fin.seek(start_pos)
		
		# Build everything we need
		self.num_chains = num_chains
		self.num_atoms = num_atoms
		self.chain = [FFEA_pdb_chain(num_atoms = num_atoms[i]) for i in range(num_chains)]

		# Read atom structures
		for c in range(self.num_chains):
			
			# Get to "ATOM" lines
			line = fin.readline()
			while True:
				if line[0:4] == "ATOM":
					fin.seek(-len(line),1)
					break
			
				line = fin.readline()

			# Read ATOM structures like a boss
			for j in range(self.chain[c].num_atoms):
				line = fin.readline()
				self.chain[c].chainID = line[21]
				self.chain[c].atom[j].set_structure(atomID=line[6:12], name=line[12:16], res=line[17:20], resID=line[22:26], occupancy=line[54:60], temperature=line[60:66], segID=line[72:76], element=line[76:78], charge=line[78:80].strip(), ffea_comment=self.get_ffea_comment(line))

		sys.stdout.write("\tdone!\n")

		#
		# Reset to beginning and read frame data
		#
		sys.stdout.write("\tReading all frame data...\n\n")
		fin.seek(start_pos)
		
		# Read all frames
		completed = False
		while not completed:
			
			# Check for EOF quickly (hopefully)
			fspos = fin.tell()
			while True:
				line = fin.readline()
				if line[0:4] == "ATOM":
					fin.seek(-len(line), 1)
					break

				if (line[0:3] == "END" and len(line) == 3) or line.strip() == "":
					completed = True
					break
			
			if completed:
				break

			# Read all data for all chains
			sys.stdout.write("\r\t\tFrames Read %d" % (self.num_frames))

			for c in range(self.num_chains):
			
				# Get to "ATOM" lines
				line = fin.readline()
				while True:

					if line[0:4] == "ATOM":
						fin.seek(-len(line),1)
						break
			
					line = fin.readline()

				# Get a frame
				frame = FFEA_frame.FFEA_frame()
				frame.num_nodes = self.num_atoms[c]
				frame.pos = [[0.0,0.0,0.0] for i in range(self.num_atoms[c])]
	
				# Read ATOM positions like a boss
				for j in range(self.chain[c].num_atoms):
					line = fin.readline()
					#print line
					frame.pos[j][0] = float(line[30:38])
					frame.pos[j][1] = float(line[38:46])
					frame.pos[j][2] = float(line[46:54])
				frame.pos = np.array(frame.pos)
				self.chain[c].frame.append(frame)
				self.chain[c].num_frames += 1
			self.num_frames += 1

			# Check finish condition
			if self.num_frames == num_frames_to_read:
				completed = True

		# Finish
		sys.stdout.write("\r\t\tFrames Read %d" % (self.num_frames))
		sys.stdout.write("\n\n\t...done! Read %d frames from file.\n" % (self.num_frames))
		fin.close()

	def get_ffea_comment(self, line):
		t = line.split("<")
		found_comment = False
		if len(t) == 2:
			t = t[1]
			if t.count(">"):
				t = t.split(">")[0]
				found_comment = True
		if found_comment:
			return t
		else:
			return ""


	def write_to_text(self, frames = None, frame_rate = 1):
		# Write frames
		if frames == None:
			frames = [0,self.num_frames]

		text = ""
		for i in range(frames[0], frames[1], frame_rate):
			#sys.stdout.write("\r\r%d frames written (%d%%)" % (i, (i * 100) / self.num_frames))
			sys.stdout.flush()
			text += ("MODEL     %4d\n" % (i + 1))
			for j in range(self.num_chains):
				for k in range(self.num_atoms[j]):
					#print j, self.num_chains, len(self.chain), k, self.num_atoms[j], len(self.chain[j].atom)
					a = self.chain[j].atom[k]
					text += ("%6s%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n" % ("ATOM  ", a.atomID, a.name, a.res, self.chain[j].chainID, a.resID, self.chain[j].frame[i].pos[k][0], self.chain[j].frame[i].pos[k][1], self.chain[j].frame[i].pos[k][2], a.occupancy, a.temperature, a.segID ,a.element, a.charge))

				text += ("TER\n")
			text+= ("ENDMDL\n")
		text += ("END\n")
		return text 


	def write_to_file(self, fname, frames = None, frame_rate = 1):

		print("Writing to " + fname + "...")

		# Write differently depending on format
		base, ext = os.path.splitext(fname)

		if ext != ".pdb":
			print_error()
			sys.stdout.write(ext + " is not a supported file extension at the moment. Defaulting to .pdb.\n")
		
		try:
			fout = open(fname, "w")
		except(IOError):
			raise IOError

		text = self.write_to_text(frames, frame_rate)
		
		fout.write(text)
		# sys.stdout.write("\r\r100% of frames written    \n")
		sys.stdout.write("flushing...")
		print("...done")
		fout.close()
	
	def build_from_traj(self, traj, scale = 1e10):

		# Reset	
		self.reset()
		
		# Chain labels
		plates = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
		
		# Each traj blob is one chain
		bin = 0
		for b in traj.blob:
			
			c = b[0]

			# Build chain
			chain = FFEA_pdb_chain()
			chain.chainID = plates[bin]

			# Sort atoms (make pseudo-atoms)
			chain.num_atoms = c.num_nodes
			for i in range(chain.num_atoms):
				atom = FFEA_pdb_atom()
				atom.set_structure()
				chain.atom.append(atom)

			# Sort frames
			chain.num_frames = traj.num_frames
			chain.frame = c.frame
			for f in chain.frame:
				f.pos *= scale

			self.chain.append(chain)
			self.num_chains += 1
			
			# Add data
			self.num_atoms.append(c.num_nodes)

			bin += 1

		self.num_frames = traj.num_frames
		self.valid = True
		self.empty = False

	def clear_position_data(self):
		
		for i in range(self.num_chains):
			self.chain[i].frame = []
			self.chain[i].num_frames = 0
		self.num_frames = 0
	
	def add_empty_frame(self):

		for i in range(self.num_chains):
			self.chain[i].add_empty_frame()

		self.num_frames += 1

	def translate(self, trans):

		for i in range(self.num_chains):
			for j in range(self.num_frames):
				try:
					self.chain[i].frame[j].translate(trans)
				except:
					print("Could not translate, likely due to formatting error in the translation vector provided")
					raise

	def rotate_chains_individually(self, rot):

		for i in range(self.num_chains):
			for j in range(self.num_frames):
				try:
					self.chain[i].frame[j].rotate(rot)
				except:
					print("Could not translate, likely due to formatting error in the translation vector provided")
					raise

	# rotate the full system according to rot, around origin_trans (default CM),
   #        but just frame findex (default 0)
	def rotate_full_system(self, rot, cent = None, findex = 0):
		# findex is which frame we move to pos
		if findex >= self.num_frames:
			print("Frame " + findex + " does not exist. Please specifiy a correct index")
			raise IndexError
		
		if cent == None:
			cent = np.array([0.0,0.0,0.0])
			for i in range(self.num_chains):
				cent += self.chain[i].frame[findex].calc_centroid() * self.num_atoms[i]
			cent *= 1.0 / sum(self.num_atoms)

		origin_trans = np.array([0.,0.,0.]) - cent
		self.translate(origin_trans)

		rot = np.array(rot)
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

		for i in range(self.num_chains):
			for j in range(len(self.chain[i].frame[findex].pos)):
				self.chain[i].frame[findex].pos[j] = np.dot(R, self.chain[i].frame[findex].pos[j])

		# Translate back
		self.translate(-1 * origin_trans)




	def set_pos(self, pos, findex = 0):
		# findex is which frame we move to pos
		if findex >= self.num_frames:
			print("Frame " + findex + " does not exist. Please specifiy a correct index")
			raise IndexError
		
		cent = np.array([0.0,0.0,0.0])
		for i in range(self.num_chains):
			cent += self.chain[i].frame[findex].calc_centroid() * self.num_atoms[i]
		cent *= 1.0 / sum(self.num_atoms)
		self.translate(np.array(pos) - cent)

	def reset(self):

		self.valid = False
		self.empty = True
		self.num_frames = 0
		self.num_chains = 0
		self.num_atoms = []
		self.chain = []

# For FFEA purposes, split object by chains (don't worry about pointers to groups of things. If you want that, use MDAnalysis)
class FFEA_pdb_chain:

	def __init__(self, num_atoms = 0, num_frames = 0):

		self.num_atoms = num_atoms
		self.atom = [FFEA_pdb_atom() for i in range(self.num_atoms)]
		self.num_frames = num_frames
		self.frame = [FFEA_frame.FFEA_frame() for i in range(self.num_frames)]
		
	def setID(self, ID):

		self.chainID = ID

	def add_empty_frame(self):
		f = FFEA_frame.FFEA_frame()
		f.num_nodes = self.num_atoms
		f.pos = np.array([[0.0,0.0,0.0] for i in range(f.num_nodes)])
		self.frame.append(f)
		self.num_frames += 1

	def set_atom_position(self, x, y, z, frame = -1, atom = None):
		if atom == None:
			raise IndexError
		
		self.frame[frame].pos[atom][0] = float(x)
		self.frame[frame].pos[atom][1] = float(y)
		self.frame[frame].pos[atom][2] = float(z)

	def reset(self):

		self.chainID = "A"
		self.num_atoms = 0
		self.num_frames = 0
		self.atom = []
		self.frame = []


class FFEA_pdb_atom:

	def __init__(self):

		self.reset()
	
	def set_structure(self, atomID=0, name="C", res="ARG", resID=0, occupancy=1.0, temperature=0.0, segID="A", element="C", charge="0", ffea_comment = ""):

		self.atomID = int(atomID)
		self.name = name
		self.res = res
		self.resID = int(resID)
		try:
			self.occupancy = float(occupancy)
		except:
			self.occupancy = 1.0

		try:
			self.temperature = float(temperature)
		except:
			self.temperature = 0.0
	
		self.segID = segID
		self.element = element
		self.charge = charge

		self.ffea_comment = ffea_comment
		
	def reset(self):

		self.atomID = 0
		self.name = "C"
		self.res = "ARG"
		self.resID = 0
		self.occupancy = 1.0
		self.temperature = 0.0
		self.segID = "A"
		self.element = "C"
		self.charge = "0"
		self.ffea_comment = ""

class FFEA_pdb_frame:

	def __init__(self, num_atoms = 0):

		self.pos = np.array([[0.0,0.0,0.0] for i in range(num_atoms)])
	
	def calc_centroid(self):
		num_atoms = len(self.pos)
		self.centroid = (1.0 / num_atoms) * np.sum(self.pos, axis = 0)
		return self.centroid

	def get_centroid(self):
		return self.centroid
	
	def translate(self, trans):

		trans = np.array(trans)
		try:
			self.pos += trans
		except:
			raise
	
	def rotate(self, rot):
		# Translate to origin
		origin_trans = np.array([0.0,0.0,0.0]) - self.calc_centroid()
		self.translate(origin_trans)

		rot = np.array(rot)
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

		for i in range(len(self.pos)):
			self.pos[i] = np.dot(R, self.pos[i])

		# Translate back
		self.translate(-1 * origin_trans)



	def set_step(self, step):
		self.step = step

	def reset(self):

		self.pos = []
		self.step = 0
