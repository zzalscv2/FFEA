import sys, os
import numpy as np
import FFEA_trajectory

class FFEA_pdb:

	def __init__(self, fname, num_frames_to_read = 100000):
	
		# Initialise
		self.reset()
		
		# Open file for reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "PDB file " + fname  + " not found. Returning empty pdb object."
			return
			
		# Start to read. Every model is a frame, every chain is a blob
		
		# So, use first frame to get the topology first of all
		findex = 0
		num_frames = 0
		bindex = 0
		num_blobs = 0
		num_atoms = []
		
		frameend = ""
		blobend = ""
		while True:
			pos = fin.tell()
			line = fin.readline()
			if "MODEL" in line:
			
				# Frames are arranged as models. This person wasn't lazy, yaaay
				frameend = "ENDMDL"
				break
				
			elif "ATOM" in line:
				
				# 1 frame, hopefully
				frameend = ""
				fin.seek(pos)
				break
				
		while True:
			line = fin.readline().strip()
			
			if line == "" or line == "END":
				break
				
			elif len(line) < 3:
				continue
							
			elif line == frameend:
			
				# Only reading first frame for now
				break
					
			elif line[0:3] == "TER":
				blobend = "TER"
				bindex += 1
				
			if line[0:4] == "ATOM":

				if bindex == num_blobs:
					num_blobs += 1
					self.blob.append(FFEA_pdb_blob())
					num_atoms.append(0)
					
				# Initialise the atom
				try:
					self.blob[bindex].atom.append(FFEA_pdb_atom())
				except(IndexError):
					break
					
				num_atoms[-1] += 1
				self.blob[bindex].num_atoms += 1
				atom_index = int(line[6:11].strip())

				# Get structure data
				atom_type = line[12:16].strip()
				res_type = line[17:20].strip()

				# If this is an ion or solvent, get rid of it
				if res_type == "SOL" or res_type == "NA" or res_type == "CL"or res_type == "K":
					self.blob[bindex].atom.pop()
					num_atoms[-1] -= 1
					self.blob[bindex].num_atoms -= 1
					continue

				chain = line[21]
				res_index = int(line[22:26].strip())

				try:
					occupancy = float(line[54:60].strip())
				except:
					occupancy = 1.00

				try:
					temp_factor = float(line[60:66].strip())
				except:
					temp_factor = 0.00

				try:
					segment = line[72:76]
				except:
					segment = "    "

				try:
					elem = line[76:78]
				except:
					elem = "C "
			
				try:
					charge = line[78:80].remove("\n")
				except:
					charge = " "

				# Ignoring position data

				# Set structure data
				self.blob[bindex].atom[num_atoms[-1] - 1].set_structure(atom_index = atom_index, atom_type = atom_type, res_type = res_type, chain = chain, res_index = res_index, occupancy = occupancy, temp_factor = temp_factor, segment = segment, elem = elem, charge = charge)
	
		self.num_blobs = num_blobs
		self.num_atoms = num_atoms
		
		# Now, load the positions only
		fin.seek(0,0)
		finished = False
		while self.num_frames < num_frames_to_read:
		
			while True:
				pos = fin.tell()
				line = fin.readline()
				
				if line.strip() == "":
					finished = True
					break
					
				elif "MODEL" in line:
			
					# Frames are arranged as models. This person wasn't lazy, yaaay
					break
				
				elif "ATOM" in line:
				
					# 1 frame, hopefully
					fin.seek(pos)
					break
			
			if finished:
				break
					
			for i in range(self.num_blobs):

				self.blob[i].add_empty_frame()
				for j in range(self.num_atoms[i]):

					line = fin.readline()
					if "ATOM" not in line:
						continue
						
					res_type = line[17:20].strip()
					
					if res_type == "SOL" or res_type == "NA" or res_type == "CL"or res_type == "K":
						continue
					else:
					
						# Get position data
						x = float(line[30:38].strip())
						y = float(line[38:46].strip())
						z = float(line[46:54].strip())
	
						self.blob[i].frame[-1].set_atomic_pos(j, x, y, z)
				
						
				if blobend == "TER":
					while fin.readline().split()[0].strip() != blobend:
						continue
				
			if frameend == "ENDMDL":
				while fin.readline().strip() != frameend:
					continue
					
			self.num_frames += 1
			sys.stdout.write("\rRead %d frames" % (self.num_frames))
			sys.stdout.flush()
			
		# Last thing, convert all to numpy
		for i in range(self.num_blobs):
			for j in range(self.num_frames):
				self.blob[i].frame[j].pos = np.array(self.blob[i].frame[j].pos)
		
	def write_to_file(self, fname, frames = []):
	
		# Build list of frames to write
		if frames == []:
			frames = range(self.num_frames)
			
		# Open file for writing
		fout = open(fname, "w")
		modindex = -1
		for i in frames:
			modindex += 1
			sys.stdout.write("\r%d%% written" % (int(modindex * 100.0 / len(frames))))
			sys.stdout.flush()
			fout.write("MODEL      %d\n" % (modindex + 1))
			for j in range(self.num_blobs):
				
				index = 0
				chain = chr(ord("A") + j) 
				for a in self.blob[j].atom:
					fout.write("ATOM  ")
					fout.write(str(a.atom_index).rjust(5))
					fout.write("  " + a.atom_type.ljust(4))
					fout.write(a.res_type.rjust(3))
					fout.write(" " + chain)
					fout.write(str(a.res_index).rjust(4))
					fout.write(" ")
					fout.write("   ")
					fout.write("%8.3f" % self.blob[j].frame[i].pos[index][0])
					fout.write("%8.3f" % self.blob[j].frame[i].pos[index][1])
					fout.write("%8.3f" % self.blob[j].frame[i].pos[index][2])
					fout.write("%6.2f" % a.occupancy)
					fout.write("%6.2f" % a.temp_factor)
					fout.write("      ")
					fout.write(a.segment.ljust(4))
					fout.write(a.elem.rjust(2))
					fout.write(a.charge.rjust(2))
					fout.write("  \n")
					index += 1
				#fout.write("TER                  " + str(self.blob[j].atom[-1].chain) + "\n")
				fout.write("TER                  " + str(chain) + "\n")
			fout.write("ENDMDL\n")
				
		fout.write("END")
		sys.stdout.write("\r100%% written\n")
		fout.close()

	def clear_position_data(self):

		self.num_frames = 0
		for b in self.blob:
			b.frame = []

	def add_frames(self, num_frames):
		
		self.num_frames = num_frames
		for b in self.blob:
			b.frame = [FFEA_pdb_frame() for i in range(self.num_frames)]

	def build_from_traj(self, traj, scale = 1.0):

		# Reset entire trajectory
		self.reset()
		
		# Check for consistency first
		for i in range(traj.num_blobs):
			if traj.num_conformations[i] != 1:
				print("Error. Cannot build a pdb from a kinetic FFEA system. Sorry.")
				return
		
		# Set properties
		self.num_blobs = traj.num_blobs
		self.num_frames = traj.num_frames
		self.blob = [FFEA_pdb_blob() for i in range(self.num_blobs)]
		for b in self.blob:
			blob_index = self.blob.index(b)
			b.num_frames = traj.num_frames
			b.num_atoms = traj.num_nodes[blob_index][0]

		# Fill up the blobs

		# Build a pseudo-structure
		for i in range(self.num_blobs):
			res_index = -1
			if i % 10 == 0:
				res_index += 1
			for j in range(self.blob[i].num_atoms):
				anatom = FFEA_pdb_atom()
				anatom.set_structure(atom_index = j, atom_type = "C", res_type = "ARG", chain = "A", res_index = res_index, occupancy = 0, temp_factor = 1.0, segment = "", elem = "C", charge = "1.0")
				self.blob[i].atom.append(anatom)

		for i in range(self.num_frames):
			for j in range(self.num_blobs):
				aframe = FFEA_pdb_frame()
				aframe.pos = traj.blob[j][0].frame[i].pos * scale
				self.blob[j].frame.append(aframe)
	
	def translate(self, shift):
		
		for i in range(len(self.blob)):
			for j in range(len(self.blob[i].frame)):
				self.blob[i].frame[j].translate(shift)
				
	def set_pos(self, pos):
	
		# Set first frame to this pos. Then, translate all other vectors by the same amount
		c = np.array([0.0,0.0,0.0])
		num_nodes = 0
		for i in range(len(self.blob)):
			c +=  len(self.blob[i].frame[0].pos) * self.blob[i].frame[0].get_centroid()
			num_nodes += len(self.blob[i].frame[0].pos)
			
		c *= 1.0 / num_nodes
		
		shift = np.array(pos) - c
		for i in range(len(self.blob)):
			for j in range(len(self.blob[i].frame)):
				self.blob[i].frame[j].translate(shift)
				
	def rotate(self, rot):

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
			
		for i in range(len(self.blob)):
			for j in range(len(self.blob[i].frame)):
			
				# Translate to origin
				origin_trans = np.array([0.0,0.0,0.0]) - self.blob[i].frame[j].get_centroid()
				self.blob[i].frame[j].translate(origin_trans)
		
				for k in range(len(self.blob[i].frame[j].pos)):
					self.blob[i].frame[j].pos[k] = np.dot(R, self.blob[i].frame[j].pos[k])
			
				# Translate back
				self.blob[i].frame[j].translate(-1 * origin_trans)
				
	def reorder_atoms(self, start = 1):

		# For each blob in the pdb
		for b in self.blob:
			end = start + b.num_atoms
			for i in range(start, end, 1):
				b.atom[i - start].atom_index = i

	def reorder_residues(self, start = 1):

		# For each blob in the pdb
		for b in self.blob:

			res_index = start		
			current_res = b.atom[0].res_index

			for i in range(b.num_atoms):

				# Are we changing residue?
				if b.atom[i].res_index != current_res:
					current_res = b.atom[i].res_index
					res_index += 1

				b.atom[i].res_index = res_index

	def add_empty_blob(self):
		self.blob.append(FFEA_pdb_blob())
		self.num_blobs += 1
		
	def reset(self):
		
		self.num_blobs = 0
		self.blob = []
		self.num_frames = 0

class FFEA_pdb_blob:

	def __init__(self):

		# Reset object
		self.reset()	

	def add_empty_frame(self):
		self.frame.append(FFEA_pdb_frame())
		self.num_frames += 1
		
	def set_pos(self, frame_index, pos):
		if self.num_atoms == len(pos):
			self.frame[frame_index].pos = pos

	def get_num_residues(self):
		
		self.num_residues = 1
		last_index = self.atom[0].res_index
		for a in self.atom[1:]:
			if a.res_index != last_index:
				last_index = a.res_index
				self.num_residues += 1

		return self.num_residues

	def reset(self):

		self.num_frames = 0
		self.frame = []
		self.atom = []
		self.num_atoms = 0
		self.num_residues = 0

class FFEA_pdb_atom:

	def __init__(self):

		# Initialise via the set_structure function
		self.set_structure()

	def set_structure(self, atom_index = 0, atom_type = "", res_type = "", chain = "", res_index = 0, occupancy = 0, temp_factor = 1.0, segment = "", elem = "", charge = ""):

		self.atom_index = atom_index
		self.atom_type = atom_type
		self.res_type = res_type
		self.chain = chain
		self.res_index = res_index
		self.occupancy = occupancy
		self.temp_factor = temp_factor
		self.segment = segment
		self.elem = elem
		self.charge = charge

class FFEA_pdb_frame:

	def __init__(self):

		self.reset()

	def reset(self):

		self.pos = []

	def get_centroid(self):
		return (1.0 / len(self.pos)) * np.sum(self.pos, axis = 0)
		
	def translate(self, trans):
		self.pos += np.array(trans)
		
	def set_pos(self, pos):
		self.translate(np.array(pos) - self.get_centroid())
	
	def set_atomic_pos(self, atom_index, x, y, z):
		
		while atom_index >= len(self.pos):
			self.pos.append([0.0,0.0,0.0])
		
		self.pos[atom_index] = [x,y,z] 
