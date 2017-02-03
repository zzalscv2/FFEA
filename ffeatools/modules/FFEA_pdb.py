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

class FFEA_pdb:

	def __init__(self, fname = "", num_frames = 100000):
	
		self.reset()

		# Empty fname give an empty object
		if fname == "":
			return

		try:
			self.load(fname, num_frames = num_frames)

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

	def load(self, fname, num_frames = 100000):

		sys.stdout.write("Loading FFEA pdb file...\n")

		# File format?
		base, ext = os.path.splitext(fname)
		if (ext == ".pdb"):
			try:
				self.load_pdb(fname, num_frames_to_read = num_frames)
			except:
				raise	
		else:
			raise FFEAIOError(fname=fname, fext=[".pdb"]) 
		
		self.valid = True
		sys.stdout.write("...done!\n")

	def load_pdb(self, fname, num_frames_to_read = 100000):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise IOError("File '" + fname + "' not found.")

				
		# Make robust against dodgy formating by the user
		sys.stdout.write("\tReading first frame to determine structure sizes...")
		models_exist = False
		chains_exist = True
		line = fin.readline()
		while(True):
			try:
				if line [0:5] == "MODEL":

					# PDB file set up for frames
					models_exist = True

					# Seek back a line
					fin.seek(-len(line),1)
					break

				elif line[0:4] == "ATOM":
					if line[21].strip() == "":
						chains_exist = False

					# Seek back a line
					fin.seek(-len(line),1)

					break

				line = fin.readline()

			except(IndexError):

				# EOF check
				if not line:
					sys.exit("EOF found before MODEL or ATOM. File badly formatted (we think...)")

		sys.stdout.write("done!\n")

		# Save pointer to the start
		start_pos = fin.tell()


		# We are ready to read the first frame, but we will not build objects, merely store the data for fast reading
		self.num_chains = 1
		self.num_frames = 1
		self.num_atoms = [0]

		if models_exist:
			fin.readline()

		ters_exist = False
		endmdls_exist = False

		sys.stdout.write("\tBuilding container objects...")
		chainIDs = []
		chainID = "A"
		while(True):

			line = fin.readline()
			if line == "":
				break
			try:
				# Different things depeinding on lines
				if line[0:4] == "ATOM":
					
					if chains_exist:
						chainID = line[21]		

					#if chainID != last_chain:
					#	if chainID in chainIDs:
					#		# No models, no ters, badly formatted. But we can deal with it
					#		break

						
					self.num_atoms[-1] += 1

				elif line[0:3] == "TER":

					chainIDs.append(chainID)
					chainID = chr(ord(chainID) + 1)					
					self.num_chains += 1
					self.num_atoms.append(0)

					# Register that ters exist
					ters_exist = True
					# Ignore. Next ATOM chain ID will catch this chain ending
				
				elif line[0:6] == "ENDMDL":
					
					endmdls_exist = True
					break
				else:
					continue
			
			except(IndexError):

				# EOF check
				if not line:

					# We will assume a single frame (rather than returning an error)
					break

		if self.num_atoms[-1] == 0:
			self.num_atoms.pop()
			self.num_chains -= 1

		fin.seek(start_pos)
		sys.stdout.write("done!\n")

		# We've got everything! Now we can read without having to check things all the time


		#
		# Build empty containers
		#

		self.chain = [FFEA_pdb_chain(num_atoms = self.num_atoms[i]) for i in range(self.num_chains)]

		# Load atom structure from first frame
		sys.stdout.write("\tLoading atomic structure from first frame...")
		if models_exist:
			fin.readline()

		for i in range(self.num_chains):

			# Chain ID first (no need to store on every atom)
			chain_start_pos = fin.tell()
			line = fin.readline()
			self.chain[i].setID(line[21])
			fin.seek(chain_start_pos)

			for j in range(self.num_atoms[i]):
				line = fin.readline()
					
				# Now atom structure
				self.chain[i].atom[j].set_structure(atomID=line[6:12], name=line[12:16], res=line[17:20], resID=line[22:26], occupancy=line[54:60], temperature=line[60:66], segID=line[72:76], element=line[76:78], charge=line[78:80])
			if ters_exist:
				fin.readline()

		
		fin.seek(start_pos)
		sys.stdout.write("done!\n")

		#
		# Load potential trajectory data
		#
		self.num_frames = 0
		lindex = 1
		sys.stdout.write("\tLoading trajectory data (if it exists)...\n")

		if num_frames_to_read == 0:
			sys.stdout.write("\t\tRequested to load 0 frames. Will upgrade to a single frame\n\n")
			num_frames_to_read = 1

		while(True):

			# Check for completion
			if self.num_frames >= num_frames_to_read:
				break

			sys.stdout.write("\r\t\tFrames read = %d" % (self.num_frames))
			sys.stdout.flush()
			try:
				if models_exist:

					# MODELS line
					line = fin.readline()
					lindex += 1
					if line.strip() == "" or line[0:3] == "END":
						break

				# Check for EOF
				frame_start_pos = fin.tell()
				line = fin.readline()
				if line.strip() == "" or line[0:3] == "END":
					break
				fin.seek(frame_start_pos)

				for i in range(self.num_chains):
				
					# Add new frame to object
					self.chain[i].add_empty_frame()

					for j in range(self.num_atoms[i]):
						line = fin.readline()
						lindex += 1

						# Get trajectory data and add to object
						try:
							self.chain[i].set_atom_position(line[30:38], line[38:46], line[46:54], atom = j)
						except(IndexError):
							print "Unable to set position for atom " + str(j) + " in chain " + str(i) + " for frame "
							raise
	
					if ters_exist:
						fin.readline()
						lindex += 1
				if endmdls_exist:
					fin.readline()
					lindex += 1

				self.num_frames += 1

			except(IndexError):

				if not line:

					# EOF
					break
				else:
					"Couldn't read file at line " + str(lindex)
					raise FFEAFormatError
		fin.close()
		sys.stdout.write("\r\t\tFrames read = %d\n\n" % (self.num_frames))
		sys.stdout.write("\t...done!\n")

	def write_to_file(self, fname):

		print "Writing to " + fname + "..."

		# Write differently depending on format
		base, ext = os.path.splitext(fname)

		if ext != ".pdb":
			print_error()
			sys.stdout.write(ext + " is not a supported file extension at the moment. Defaulting to .pdb.\n")
		
		try:
			fout = open(fname, "w")
		except(IOError):
			"Cannot open " + fname + " for writing."
			raise
		
		for i in range(self.num_frames):
			sys.stdout.write("\r\r%d frames written (%d%%)" % (i, (i * 100) / self.num_frames))
			sys.stdout.flush()
			fout.write("MODEL     %4d\n" % (i + 1))
			for j in range(self.num_chains):
				for k in range(self.num_atoms[j]):
					a = self.chain[j].atom[k]
					fout.write("%6s%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n" % ("ATOM  ", a.atomID, a.name, a.res, self.chain[j].chainID, a.resID, self.chain[j].frame[i].pos[k][0], self.chain[j].frame[i].pos[k][1], self.chain[j].frame[i].pos[k][2], a.occupancy, a.temperature, a.segID ,a.element, a.charge))

				fout.write("TER\n")
			fout.write("ENDMDL\n")
		fout.write("END\n")
		sys.stdout.write("\r\r100%% of frames written\n")
		print "...done"
		fout.close()
	
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
		self.frame = [FFEA_pdb_frame() for i in range(self.num_frames)]
		
	def setID(self, ID):

		self.chainID = ID

	def add_empty_frame(self):
		
		self.frame.append(FFEA_pdb_frame(num_atoms = self.num_atoms))
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
	
	def set_structure(self, atomID=0, name="C", res="ARG", resID=0, occupancy=1.0, temperature=0.0, segID="A", element="C", charge="0"):

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
		
	def reset(self):

		self.pos = []
