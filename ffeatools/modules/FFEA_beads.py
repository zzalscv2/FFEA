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
from FFEA_exceptions import *
import FFEA_pdb

class FFEA_beads:
 
	def __init__(self, fname = "", motion_state = "STATIC", topology = None, node = None):

		self.reset()
		if fname == "":
			self.valid = True
			sys.stdout.write("Empty beads object initialised\n")
			return

		#try:
		self.load(fname, motion_state, topology, node)
		#except:
		#	raise

	def load(self, fname, motion_state, topology, node):
		sys.stdout.write("Loading FFEA beads...\n")

		# load the PDB
		base, ext = os.path.splitext(fname)
		try:
			if (ext == ".pdb"):
				self.pdb = FFEA_pdb.FFEA_pdb(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".pdb"])

		except:
			raise

		# now assign beads to elements... but we only have 
		#  elements if motion_state == DYNAMIC
		if motion_state == "DYNAMIC":
			self.assign_beads(topology, node)


		self.valid = True
		self.empty = False
		sys.stdout.write("...done!\n")


	def get_beads_preassignment(self, text):
		beads_assignment = []
		l = text.split("=")
		if len(l) != 2: return beads_assignment
		if l[0].strip() != "nodes": return beads_assignment
		else:
			for i in l[1].split(","):
				beads_assignment.append(int(i))
		return beads_assignment
			

	def assign_beads(self, topology, node):

		# initialise b_elems with -1
		for i in range(self.pdb.num_chains):
			aux = []
			for j in range(self.pdb.chain[i].num_atoms):
				aux.append(-1)
			self.b_elems.append(aux)
				

		for i in range(self.pdb.num_chains):
			for j, a in enumerate(self.pdb.chain[i].atom):
				d2_0 = 1e9
				beads_assignment = self.get_beads_preassignment(a.ffea_comment)
				# get the closest element to this bead:
				for e_ndx, e in enumerate(topology.element):
					work = True
					if len(beads_assignment) > 0:
						work = False
						for ba in beads_assignment:
							if e.n.count(ba):
								work = True
								break
					if (work == False): continue
					ecm = e.calc_centroid(node)
					d2 = (ecm[0] - self.pdb.chain[i].frame[0].pos[j][0])**2 + \
					     (ecm[1] - self.pdb.chain[i].frame[0].pos[j][1])**2 + \
					     (ecm[2] - self.pdb.chain[i].frame[0].pos[j][2])**2
							
					if (d2 < d2_0):
						d2_0 = d2
						self.b_elems[i][j] = e_ndx

				#print a.name, a.res, self.pdb.chain[i].frame[0].pos[j][0:3],\
                  #a.ffea_comment, self.b_elems[i][j]
				if self.b_elems[i][j] == -1: 
					raise

	def reset(self):
		self.pdb = None
		self.b_elems = [] # element indices, double list for [chain][atom].


		
