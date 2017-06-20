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

import numpy as np
import sys, os
import FFEA_pdb, FFEA_node, FFEA_trajectory

class FFEA_kinetic_map:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print("Error. Map File " + fname  + " not found.")
			return

		# Header
		# Do we need to convert the file?
		line = fin.readline()
		if "Dense" in line:
	
			# Convert map to sparse first
			scriptdir = os.path.dirname(os.path.realpath(__file__))
			print("Converting dense map to sparse...")
			base, ext = os.path.splitext(os.path.abspath(fname))
			os.system("python " + scriptdir + "/../../FFEA_initialise/FFEA_mapping_tools/Kinetic_FFEA_map_convert_to_sparse.py " + fname + " " + base + "_sparse" + ext)
			fname = base + "_sparse" + ext
			print("done!")

		elif "Sparse" not in line:
			sys.exit("Error. This may not be an FFEA kinetic map. Expected '(Dense)' or '(Sparse)' in first line of file.")

		fin.close()
		fin = open(fname, "r")
		if fin.readline().strip() != "FFEA Kinetic Conformation Mapping File (Sparse)":
			self.reset()
			print("Error. Incorret header layout.")
			return
 
		self.num_columns = int(fin.readline().split()[1])
		self.num_rows = int(fin.readline().split()[1])
		self.num_entries = int(fin.readline().split()[1])

		self.entry = [0.0 for i in range(self.num_entries)]
		self.key = [0 for i in range(self.num_rows + 1)]
		self.col = [0 for i in range(self.num_entries)]

		if fin.readline().strip() != "map:":
			self.reset()
			print("Error. Incorret header layout.")
			return

		# Read entries
		print("Reading entries...")
		i = -1
		for entry in fin.readline().split()[2:]:
			i += 1
			self.entry[i] = float(entry)
		print("done!")

		# Read key
		print("Reading key...")
		i = -1
		for entry in fin.readline().split()[2:]:
			i += 1
			self.key[i] = int(entry)
		print("done!")

		# Read columns
		i = -1
		print("Reading columns...")
		for entry in fin.readline().split()[2:]:
			i += 1
			self.col[i] = int(entry)
		print("done!")
		print("Map reading completed. Ready for application.")
		return

	def apply_sparse(self, base):
		
		# Get base type and apply appropriately
		if isinstance(base, FFEA_pdb.FFEA_pdb):
			basetype = "pdb"
			if self.num_columns != len(base.blob[0].frame[0].pos):
				print("Error. Map expects " + str(self.num_columns) + " nodes, but found " + str(len(base.blob[0].frame[0].pos)))
				return
			else:
				total_new_nodes = []
				count = 0
				for f in base.blob[0].frame:
					count += 1
					new_nodes = np.array([[0.0,0.0,0.0] for i in range(self.num_rows)])
					for i in range(self.num_rows):
						for j in range(self.key[i], self.key[i + 1]):
							new_nodes[i] += self.entry[j] * f.pos[self.col[j]]
							
					print(str(count) + " frames calculated")
					total_new_nodes.append(new_nodes)
						
			return total_new_nodes

		elif isinstance(base, FFEA_node.FFEA_node) or isinstance(base, FFEA_frame.FFEA_frame):
			basetype = "FFEA_node"
			if self.num_columns != len(base.pos):
				print("Error. Map expects " + str(self.num_columns) + " nodes, but found " + str(len(base.pos)))
				return
			else:
				new_nodes = np.array([[0.0,0.0,0.0] for i in range(self.num_rows)])
				for i in range(self.num_rows):
					for j in range(self.key[i], self.key[i + 1]):
						new_nodes[i] += self.entry[j] * base.pos[self.col[j]]
			return new_nodes

		elif isinstance(base, FFEA_trajectory.FFEA_traj_blob_frame):
			basetype = "FFEA_traj_blob_frame"
			if self.num_columns != len(base.pos):
				print("Error. Map expects " + str(self.num_columns) + " nodes, but found " + str(len(base.pos)))
				return
			else:
				new_nodes = np.array([[0.0,0.0,0.0] for i in range(self.num_rows)])
				for i in range(self.num_rows):
					for j in range(self.key[i], self.key[i + 1]):
						new_nodes[i] += self.entry[j] * base.pos[self.col[j]]
			
		return new_nodes

	def reset(self):
		self.num_entries = 0
		self.num_rows = 0
		self.num_columns = 0
		self.entry = []
		self.key = []
		self.col = []
