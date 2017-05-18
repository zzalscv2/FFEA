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

class FFEA_pin:

	def __init__(self, fname = ""):
	
		self.reset()

		# Empty fname give an empty object
		if fname == "":
			return

		try:
			self.load(fname)

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

	def load(self, fname):

		sys.stdout.write("Loading FFEA pin file...")

		# File format?
		base, ext = os.path.splitext(fname)
		try:
			if ext == ".pin":
				self.load_pin(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".pin"])

		except:
			raise

		self.valid = True
		sys.stdout.write("done!\n")

	def load_pin(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea pinned nodes file" and line != "walrus pinned nodes file":
			print("\tExpected 'ffea pinned nodes file' but found " + line)
			raise TypeError

		num_pinned_nodes = int(fin.readline().split()[1])

		fin.readline()

		# Read pinned nodes now
		pintype = 0	
		while(True):
			line = fin.readline().strip()
			if line == "":
				break
			else:
				self.add_pinned_node(line)

		fin.close()

	def add_pinned_node(self, anint):

		self.index.append(int(anint))
		self.num_pinned_nodes += 1
		
	def remove_pinned_node(self, anint):
		
		try:
			self.index.remove(anint)
			self.num_pinned_nodes -= 1

		except(ValueError):
			print("Index " + str(anint) + " not in list.")
			

	def print_details(self):

		print "num_pinned_nodes = %d" % (self.num_pinned_nodes)
		sleep(1)

		outline = ""
		for i in self.index:
			outline += str(i) + " "
			
		print outline
	
	def write_to_file(self, fname):
		
		with open(fname, "w") as f:
			f.write("ffea pinned nodes file\nnum_pinned_nodes %d\npinned nodes:\n" % (self.num_pinned_nodes))
			for i in self.index:
				f.write("%d\n" % (i))


	def pin_radially(self, node, oindex, radius, top=None, linear=0):
		
		# Reset first
		self.reset()

		# Pin all within radius
		origin = node.pos[oindex]
		
		nindex = -1

		# Get relevent index list
		if linear == 0:
			indices = range(node.num_nodes)
		else:
			if top == None:
				print "Linear indices cannot be found without a topology. Defaulting to all nodes..."
				range(node.num_nodes)
			else:
				indices = []
				for el in top.element:
					for i in el.n[0:4]:
						indices.append(i)

				indices = list(set(indices))
				
		for i in indices:
			
			d = np.linalg.norm(node.pos[i] - origin)
			if d < radius:
				self.add_pinned_node(i)

			
	def reset(self):

		self.index = []
		self.num_pinned_nodes = 0
		self.valid = False
