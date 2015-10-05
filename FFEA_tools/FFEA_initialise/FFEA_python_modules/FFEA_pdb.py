import sys, os
import numpy as np

class FFEA_pdb:

	def __init__(self, fname):

		# Initialise variables
		self.atom = []
		self.num_atoms = 0

		# Open and read file
		try:
			fin = open(fname, "r")
		except(IOError):
			return

		lines = fin.readlines()
		fin.close()

		for line in lines:

			# Only get a single frame for now (very basic, update me continuously)
			if line[0:6] == "ENDMDL":
				break

			if line[0:4] == "ATOM":
				self.num_atoms += 1
				atom_index = self.num_atoms
				#atom_index = int(line[6:11].strip())
				atom_type = line[12:16].strip()
				residue_type = line[17:20].strip()
				chain = line[21]
				residue_index = int(line[22:26].strip())
				x = float(line[30:38].strip())
				y = float(line[38:46].strip())
				z = float(line[46:54].strip())
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

				self.atom.append(FFEA_pdb_atom(atom_index, atom_type, residue_type, chain, residue_index, np.array([x,y,z]), occupancy, temp_factor, segment, elem, charge))

		self.calc_centroid()

	def write_to_file(self, fname):

		# Open file
		fout = open(fname, "w")	
		for a in self.atom:
			fout.write("ATOM  ")
			fout.write(str(a.atom_index).rjust(5))
			fout.write("  " + a.atom_type.ljust(4))
			fout.write(a.residue_type.rjust(3))
			fout.write(" " + a.chain)
			fout.write(str(a.residue_index).rjust(4))
			fout.write(" ")
			fout.write("   ")
			fout.write("%8.3f" % a.pos[0])
			fout.write("%8.3f" % a.pos[1])
			fout.write("%8.3f" % a.pos[2])
			fout.write("%6.2f" % a.occupancy)
			fout.write("%6.2f" % a.temp_factor)
			fout.write("      ")
			fout.write(a.segment_id.ljust(4))
			fout.write(a.element.rjust(2))
			fout.write(a.charge.rjust(2))
			fout.write("  \n")
		fout.write("TER                  " + str(self.atom[-1].chain))
		fout.close()

	def scale(self, scale):

		for a in self.atom:
			for j in range(3):
				a.pos[j] *= scale
			
	# Removes and returns chuck of pdb. Keeps the leftover
	def extract_atoms(self, atomlist):
		
		atomlist.sort()

		fromatom = atomlist[0]
		toatom = atomlist[-1]

		if fromatom > toatom:
			print "Lower index > upper index. Switching..."
			fromatom, toatom = toatom, fromatom

		# Check for whole residues
		print "Specified extraction of atom " + str(fromatom) + " to " + str(toatom) + "..."
		while(True):
			try:
				if self.atom[fromatom].residue_index != self.atom[fromatom - 1].residue_index:
					break
				else:
					fromatom -= 1

			except(IndexError):
				break

		while(True):
			try:
				if self.atom[toatom].residue_index != self.atom[toatom + 1].residue_index:
					break

				toatom += 1

			except(IndexError):
				break

		print "Checked for complete residues. This corresponds to " + str(self.atom[fromatom].atom_index) + " up to and including " + str(self.atom[toatom - 1].atom_index)

		# Creating the new object
		chunkpdb = FFEA_pdb("")
		
		# Appending to new, removing from old
		for i in atomlist:
			chunkpdb.add_atom(self.atom[i])

		self.remove_atoms(atomlist)
		
		return chunkpdb

	def get_atomlist_radially(self, origin, radius):
		
		atomlist = []
		for i in range(len(self.atom)):
			distance = np.linalg.norm(self.atom[i].pos - origin)
			if distance < radius:
				atomlist.append(i)

		return atomlist
	
	def add_atom(self, ffeaatom):
		self.atom.append(ffeaatom)
		self.num_atoms += 1
		self.atom[-1].atom_index = self.num_atoms

	def remove_atoms(self, atomlist):
		offset = len(atomlist)
		indexlist = range(len(self.atom))
		indexlist.reverse()
		for i in indexlist:
			if i in atomlist:
				offset -= 1
				self.atom.pop(i)
			else:
				self.atom[i].atom_index -= offset

		self.num_atoms -= len(atomlist)

	def add_pdb(self, ffeapdb):
		for a in ffeapdb.atom:
			self.atom.append(a)
			self.num_atoms += 1
			self.atom[-1].atom_index = self.num_atoms


	def calc_centroid(self):
		self.centroid = np.array([0.0, 0.0, 0.0])
		for a in self.atom:
			self.centroid += a.pos

		self.centroid *= 1.0/self.num_atoms

	def translate(self, translation):
		for a in self.atom:
			a.pos += translation

	def rotate(self, point, axis, angle):
		
		# First make sure axis is normalised!
		mag = np.linalg.norm(axis)
		axis *= 1.0/mag

		# And angle is radians
		angle = np.radians(angle)

		# Move point to origin
		self.translate(-1 * point)
		
		# Now do rotation
		c = np.cos
		s = np.sin
		r00 = c(angle) + axis[0]**2 * (1 - c(angle))
		r01 = axis[0] * axis[1] * (1 - c(angle)) - axis[2] * s(angle)
		r02 = axis[0] * axis[2] * (1 - c(angle)) + axis[1] * s(angle)
		r10 = axis[0] * axis[1] * (1 - c(angle)) + axis[2] * s(angle)
		r11 = c(angle) + axis[1]**2 * (1 - c(angle))
		r12 = axis[1] * axis[2] * (1 - c(angle)) - axis[0] * s(angle)
		r20 = axis[0] * axis[2] * (1 - c(angle)) - axis[1] * s(angle)
		r21 = axis[1] * axis[2] * (1 - c(angle)) + axis[0] * s(angle)
		r22 = c(angle) + axis[2]**2 * (1 - c(angle))
		R = np.array([[r00,r01,r02],[r10,r11,r12],[r20,r21,r22]])
		
		for a in self.atom:
			a.pos = np.dot(R, a.pos)
			
		# Move back
		self.translate(point)


class FFEA_pdb_atom:

	def __init__(self, ai, at, rt, c, ri, p, oc, temp, seg, elem, charge):
		self.set_params(ai, at, rt, c, ri, p, oc, temp, seg, elem, charge)

	def set_params(self, ai, at, rt, c, ri, p, oc, temp, seg, elem, charge):
		self.atom_index = ai
		self.atom_type = at
		self.residue_index = ri
		self.residue_type = rt
		self.chain = c
		self.pos = p
		self.occupancy = oc
		self.temp_factor = temp
		self.segment_id = seg
		self.element = elem
		self.charge = charge

def write_frames_to_file(pdb, pdbfname):

	# Open file
	fout = open(pdbfname, "w")
	
	# Write out each frame as a model
	for frame in pdb:
		fout.write("MODEL      %d\n" % (pdb.index(frame) + 1))
		for a in frame.atom:
			fout.write("ATOM  ")
			fout.write(str(a.atom_index).rjust(5))
			fout.write("  " + a.atom_type.ljust(4))
			fout.write(a.residue_type.rjust(3))
			fout.write(" " + a.chain)
			fout.write(str(a.residue_index).rjust(4))
			fout.write(" ")
			fout.write("   ")
			fout.write("%8.3f" % a.pos[0])
			fout.write("%8.3f" % a.pos[1])
			fout.write("%8.3f" % a.pos[2])
			fout.write("%6.2f" % a.occupancy)
			fout.write("%6.2f" % a.temp_factor)
			fout.write("      ")
			fout.write(a.segment_id.ljust(4))
			fout.write(a.element.rjust(2))
			fout.write(a.charge.rjust(2))
			fout.write("  \n")
		fout.write("TER\nENDMDL\n")
	fout.close()
