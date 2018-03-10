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
import numpy as np
from FFEA_exceptions import *
import FFEA_pdb

class FFEA_beads:
 
	def __init__(self, fname = "", motion_state = "STATIC", scale = 1.0, topology = None, node = None, assignBeads = True):

		self.reset()
		self.scale = scale
		if fname == "":
			self.valid = True
			sys.stdout.write("Empty beads object initialised\n")
			return

		#try:
		self.load(fname, motion_state, topology, node, assignBeads)
		#except:
		#	raise

	def load(self, fname, motion_state, topology, node, assignBeads):
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

		# and rescale the positions:
		self.rescale(self.scale)

		# now assign beads to elements... but we only have 
		#  elements if motion_state == DYNAMIC
		if motion_state == "DYNAMIC" and assignBeads:
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
			

	# initialise_assignment should have been called
	def assign_beads(self, topology, node):

		# initialise b_elems with -1
		for i in range(self.pdb.num_chains):
			aux = [-1] * self.pdb.chain[i].num_atoms
			self.b_elems.append(aux)
				
		# precalculate the CMs of every element.
		eCM = []
		for e_ndx, e in enumerate(topology.element):
			eCM.append( e.calc_centroid(node) ) 
		eCM = np.array(eCM)

		# put all the atoms into a numpy array:
		npA = []
		for i in range(self.pdb.num_chains): # new 
			for j, a in enumerate(self.pdb.chain[i].atom):
				npA.append(self.pdb.chain[i].frame[0].pos[j])
		npA = np.array(npA) 
		# xyD2[e,a] is the squared distance matrix 
      #           between atoms and element centroids
		xyD2 = ((eCM[:,:,None] - npA[:,:,None].T) ** 2 ).sum(1)

		a_cnt = -1 # new 
		for i in range(self.pdb.num_chains):
			for j, a in enumerate(self.pdb.chain[i].atom):
				a_cnt += 1 # new 
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
					d2 = xyD2[e_ndx, a_cnt]
							
					if (d2 < d2_0):
						d2_0 = d2
						self.b_elems[i][j] = e_ndx

				if self.b_elems[i][j] == -1: 
					print "ABORTING: an element could not be assigned"
					raise


	def rescale(self, factor):
		for i in range(self.pdb.num_chains):
			for j, a in enumerate(self.pdb.chain[i].atom):
				for d in range(3):
					self.pdb.chain[i].frame[0].pos[j][d] *= factor


	def reset(self):
		self.valid = False
		self.empty = True
		self.pdb = None
		self.b_elems = [] # element indices, double list for [chain][atom].
		self.scale = 1.0


		
