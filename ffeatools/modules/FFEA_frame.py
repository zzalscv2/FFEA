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
import FFEA_node

class FFEA_frame(FFEA_node.FFEA_node):

	# Load with file object already open
	def load_from_traj(self, fo):
		
		start = fo.tell()
		
		while(True):
			prev = fo.tell()
			line = fo.readline().split()
			try:
				self.pos.append([float(line[i]) for i in range(3)])
				self.vel.append([float(line[i]) for i in range(3,6)])
			except(IndexError):
			
				# EOF
				fo.seek(start)
				return 1

			except(ValueError):
				if line[0] == "*" or line[0] == "Blob":
				
					# Ready for next frame
					fo.seek(prev)
					break
				else:
					# Halfway through a written frame
					fo.seek(start)
					return 1
					
		# Numpy it up for speed
		self.pos = np.array(self.pos)
		self.vel = np.array(self.vel)
		self.num_nodes = len(self.pos)
		self.num_surface_nodes = self.num_nodes
		
		return 0
		
	# Load only positions, knowing the number of nodes, 
   #      from file object already open
	def load_from_traj_faster(self, fo):

		if self.num_nodes == 0:
			self.load_from_traj(fo)
			return 0
		
		start = fo.tell()
		
		self.pos = np.empty([self.num_nodes, 3])
		self.vel = np.empty([self.num_nodes, 3])
		cnt = 0 
		while(True):
			prev = fo.tell()
			line = fo.readline().split()
			try:
				for i in range(3):
					self.pos[cnt,i] = float(line[i])
					self.vel[cnt,i] = float(line[i+3])
				cnt += 1
			except(IndexError):
			
				# EOF
				fo.seek(start)
				return 1

			except(ValueError):
				if line[0] == "*" or line[0] == "Blob":
				
					# Ready for next frame
					fo.seek(prev)
					break
				else:
					# Halfway through a written frame
					fo.seek(start)
					return 1
					
		self.num_surface_nodes = self.num_nodes
		return 0
		
	# Load only positions, knowing the number of nodes, 
   #      from file object already open
	def load_from_traj_onlynodes_faster(self, fo):

		if self.num_nodes == 0:
			self.load_from_traj(fo)
			return 0
		
		start = fo.tell()
		
		self.pos = np.empty([self.num_nodes, 3])
		cnt = 0 
		try:
		   #read fixed sized lines
		   # text = fo.read(line_length*self.num_nodes).split("\n")
		   for ln in range(self.num_nodes):
		      self.pos[ln] = np.fromstring(fo.readline(), dtype=float, sep=" ", count=3)
		except:
		   fo.seek(start)
		   return 1

		self.num_surface_nodes = self.num_nodes
		return 0
		
	def build_from_node(self, node):
	
		self.num_nodes = node.num_nodes
		self.num_surface_nodes = node.num_surface_nodes
		self.num_interior_nodes = node.num_interior_nodes
		self.pos = node.pos
		
	def write_to_traj(self, fo):

		for p in self.pos:
			fo.write("%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n" % (p[0], p[1], p[2], 0, 0, 0, 0, 0, 0, 0))

	# Function to calculate normals at each node, average of connecting faces
	def calc_normals(self, surf):

		self.normal = np.array([[0.0, 0.0, 0.0] for i in range(len(self.pos))])
		for f in surf.face:
			norm = f.get_normal(self)
			for n in f.n:
				self.normal[n] += norm

		for n in self.normal:
			n /= np.linalg.norm(n)

	def set_step(self, step):
		self.step = step

		
	def reset(self):
		self.num_nodes = 0
		self.num_surface_nodes = 0
		self.num_interior_nodes = 0
		self.step = 0
		self.pos = []
		self.vel = []
		self.normal = []
