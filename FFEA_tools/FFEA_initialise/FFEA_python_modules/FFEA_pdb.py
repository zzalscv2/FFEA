import sys, os
import numpy as np

class FFEA_pdb:

	def __init__(self, fname):

		# Initialise everything
		self.reset()

		# Open file for reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. File " + fname  + " not found."
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
			if line[0:3] == "END":
				break

			# Not interested in this
			elif len(line.strip()) < 5:
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
				#atom_index = self.num_atoms
				atom_index = int(line[6:11].strip())

				# Get structure data
				atom_type = line[12:16].strip()
				res_type = line[17:20].strip()
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
		fout.close()
	
	def clear_position_data(self):

		self.num_frames = 0
		for b in self.blob:
			b.frame = []

	def add_frames(self, num_frames):
		
		self.num_frames = num_frames
		for b in self.blob:
			b.frame = [FFEA_pdb_frame() for i in range(self.num_frames)]

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

	def reset(self):

		self.num_frames = 0
		self.frame = []
		self.atom = []
		self.num_atoms = 0

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


	def set_atomic_pos(self, atom_index, x, y, z):
		
		while atom_index >= len(self.pos):
			self.pos.append([0.0,0.0,0.0])
		
		self.pos[atom_index] = [x,y,z] 