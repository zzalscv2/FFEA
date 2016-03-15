import sys, os
import numpy as np
import FFEA_trajectory

class FFEA_pdb:

	def __init__(self, fname, num_frames_to_read = 1000000):

		# Initialise everything
		self.reset()

		# Open file for reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "PDB file " + fname  + " not found. Returning empty pdb object."
			return

		# The file may contain frames as models, models of separate structures, or just a single structure. 
		# Work this out first with an initial sweep
		print("Making an initial sweep of the structure...")
		num_models = 0
		num_atoms = []
		lines = fin.readlines()
		for line in lines:

			# Not interested in this
			if len(line.strip()) < 5:
				continue
				
			# Have we got models?
			if line[0:5] == "MODEL":
				num_models += 1
				num_atoms.append(0)

			elif line[0:4] == "ATOM":

				# If not in models
				if num_models == 0 and len(num_atoms) == 0:
					num_atoms.append(0)

				num_atoms[-1] += 1

		print("done!")
		print("num_models = " + str(num_models))

		# This should be fixed for multiple blobs with multiple frames in future
		if len(set(num_atoms)) == 1:

			print("num_atoms = " + str(num_atoms[0]) + " for all models. Models are frames.")
			if num_models == 0:
				self.num_frames = 1
			else:
				self.num_frames = num_models
			self.num_blobs = 1
		else:
			print("num_atoms is different for all models. Models are different structures.")
			self.num_blobs = num_models
			self.num_frames = 1

		if num_frames_to_read < self.num_frames:
			self.num_frames = num_frames_to_read

		# Setup the frames and structures
		if num_models == 0:
			self.num_blobs = 1
		self.blob = [FFEA_pdb_blob() for i in range(self.num_blobs)]

		# Load in the structure and trajectory
		print("Reading structure and trajectory...")
		num_atoms = 0
		model_index = 0
		atom_index = 0
		blob_index = 0
		frame_index = 0
	
		# Get a single frame for everything
		for b in self.blob:
			b.frame.append(FFEA_pdb_frame())

		for line in lines:


			# Finished
			if line.strip() == "END" or self.num_frames == frame_index:
				break

			# Not interested in this
			elif len(line.strip()) < 1:
				continue
			
			elif line[0:6] == "ENDMDL":
				model_index += 1

				# May need fixing for multiple blobs with multiple frames
				if self.num_frames == 1:
					blob_index = model_index
				else:
					frame_index = model_index

					# Might need another frame
					if frame_index != self.num_frames:
						self.blob[blob_index].frame.append(FFEA_pdb_frame())

			elif line[0:3] == "TER":
				self.blob[blob_index].num_atoms = num_atoms
				atom_index = 0
				num_atoms = 0

			elif line[0:4] == "ATOM":

				# Initialise the atom
				self.blob[blob_index].atom.append(FFEA_pdb_atom())
				atom_index = int(line[6:11].strip())

				# Get structure data
				atom_type = line[12:16].strip()
				res_type = line[17:20].strip()

				# If this is an ion or solvent, get rid of it
				if res_type == "SOL" or res_type == "NA" or res_type == "CL"or res_type == "K":
					self.blob[blob_index].atom.pop()
					self.blob[blob_index].num_atoms -= 1
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

				# Get position data
				x = float(line[30:38].strip())
				y = float(line[38:46].strip())
				z = float(line[46:54].strip())

				# Set structure data
				self.blob[blob_index].atom[num_atoms].set_structure(atom_index = atom_index, atom_type = atom_type, res_type = res_type, chain = chain, res_index = res_index, occupancy = occupancy, temp_factor = temp_factor, segment = segment, elem = elem, charge = charge)

				# Set position data
				self.blob[blob_index].frame[frame_index].set_atomic_pos(num_atoms, x, y, z)

				# Increment index
				num_atoms += 1

		for b in self.blob:

			# If no 'TER's
			if b.num_atoms == 0:
				b.num_atoms = len(b.atom)

			for f in b.frame:
				f.pos = np.array(f.pos)

		print("done!")
		
	def write_to_file(self, fname):

		# Open file for writing
		fout = open(fname, "w")

		# Very simple output, no headers or any stuff like that. Frames are models
		for i in range(self.num_frames):
			sys.stdout.write("\r%d%% written" % (int(i * 100.0 / self.num_frames)))
			sys.stdout.flush()
			for b in self.blob:
				fout.write("MODEL\n")
				index = 0
				for a in b.atom:
					fout.write("ATOM  ")
					fout.write(str(a.atom_index).rjust(5))
					fout.write("  " + a.atom_type.ljust(4))
					fout.write(a.res_type.rjust(3))
					fout.write(" " + a.chain)
					fout.write(str(a.res_index).rjust(4))
					fout.write(" ")
					fout.write("   ")
					fout.write("%8.3f" % b.frame[i].pos[index][0])
					fout.write("%8.3f" % b.frame[i].pos[index][1])
					fout.write("%8.3f" % b.frame[i].pos[index][2])
					fout.write("%6.2f" % a.occupancy)
					fout.write("%6.2f" % a.temp_factor)
					fout.write("      ")
					fout.write(a.segment.ljust(4))
					fout.write(a.elem.rjust(2))
					fout.write(a.charge.rjust(2))
					fout.write("  \n")
					index += 1
				fout.write("TER                  " + str(b.atom[-1].chain) + "\n")
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

	def set_single_frame(self, index):

		for b in self.blob:
			try:
				f = b.frame[index]
			except(IndexError):
				print("Frame " + str(index) + " does  not exist.\n")
				return
			
			b.frame = []
			b.frame.append(f)
			b.num_frames = 1

		self.num_frames = 1

	def scale(self, factor):

		for b in self.blob:
			for f in b.frame:
				f.pos *= factor

	def translate(self, shift):

		for b in self.blob:
			for f in b.frame:
				f.pos += np.array(shift)

	def move(self, target, frame = 0):

		c = self.calc_centroid(frame)
		self.translate(target - c)
	
	def calc_centroid(self, frame):

		c = np.array([0,0,0])
		total_num_atoms = 0

		for b in self.blob:
			c += b.frame[frame].calc_centroid() * b.num_atoms
			total_num_atoms += b.num_atoms

		c /= total_num_atoms
		return c

	def reset(self):
		
		self.num_blobs = 0
		self.blob = []
		self.num_frames = 0

class FFEA_pdb_blob:

	def __init__(self):

		# Reset object
		self.reset()	

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

	def calc_centroid(self):

		return np.mean(self.pos, axis=0)

	def set_atomic_pos(self, atom_index, x, y, z):
		
		while atom_index >= len(self.pos):
			self.pos.append([0.0,0.0,0.0])
		
		self.pos[atom_index] = [x,y,z] 
