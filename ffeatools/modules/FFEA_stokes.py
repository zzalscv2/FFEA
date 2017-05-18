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
from numpy import pi

class FFEA_stokes:

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

		sys.stdout.write("Loading FFEA stokes file...")

		# File format?
		base, ext = os.path.splitext(fname)
		try:
			if ext == ".stokes":
				self.load_stokes(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".stokes"])

		except:
			raise

		self.valid = True
		sys.stdout.write("done!\n")


	def load_stokes(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea stokes radii file" and line != "walrus stokes radii file":
			print("\tExpected 'ffea stokes radii file' but found " + line)
			raise TypeError

		num_nodes = int(fin.readline().split()[1])

		# Read stokes radii now
		while(True):
			line = fin.readline().strip()
			if line == "":
				break
			else:
				self.add_node(line)

		fin.close()

	def default(self, num_nodes, top, rad):

		self.num_nodes = num_nodes

		# Default
		self.radius = [0.0 for i in range(self.num_nodes)]

		for el in top.element:
			
			# Only linear nodes
			for n in el.n[0:4]:
				self.radius[n] = float(rad)

	def add_node(self, afloat):

		self.radius.append(float(afloat))
		self.num_nodes += 1

	def write_to_file(self, fname):

		with open(fname, "w") as f:
			f.write("ffea stokes radii file\nnum_nodes %d\n" % (self.num_nodes))
			for rad in self.radius:
				f.write("%6.3f\n" % (rad))
	
	def print_details(self):

		print "num_nodes = %d" % (self.num_nodes)
		sleep(1)

		outline = ""
		for rad in self.radius:
			outline += "%6.3f\n" % (rad)
			
		print outline
	
	def calc_drag(self, viscosity, scale = 1.0):
		
		drag = 0.0
		for r in self.radius:
			drag += 6 * pi * viscosity * r * scale

		return drag

	def reset(self):

		self.radius = []
		self.num_nodes = 0
