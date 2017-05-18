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

class FFEA_vdw:

	def __init__(self, fname = ""):
	
		self.reset()

		if fname == "":
			return

		try:
			self.load(fname)

		except FFEAFormatError as e:
			self.reset()
			print_error()
			print "Formatting error at line " + e.lin + "\nLine(s) should be formatted as follows:\n\n" + e.lstr
			raise

		except FFEAIOError as e:
			self.reset()
			print_error()
			print "Input error for file " + e.fname
			if e.fext != [""]:
				print "       Acceptable file types:"
				for ext in e.fext:
					print "       " + ext
		except IOError:
			raise

	def load(self, fname):

		sys.stdout.write("Loading FFEA vdw file...")

		# File format?
		base, ext = os.path.splitext(fname)
		try:
			if ext == ".vdw":
				self.load_vdw(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".vdw"])

		except:
			raise
	
		self.valid = True
		sys.stdout.write("done!\n")

	def load_vdw(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea vdw file" and line != "walrus vdw file":
			raise TypeError("\tExpected 'ffea vdw file' but found " + line)

		num_faces = int(fin.readline().split()[1])

		fin.readline()

		# Read vdw radii now
		while(True):
			line = fin.readline().strip()
			if line == "":
				break
			else:
				self.add_face(line)

		fin.close()

	def default(self, num_faces):
		self.num_faces = num_faces
		self.index = [-1 for i in range(num_faces)]
		
	def add_face(self, anint):

		self.index.append(int(anint))
		self.num_faces += 1
	
	def set_index(self, findex, vdwindex):
		
		try:
			self.index[findex] = int(vdwindex)
		except:
			raise
	
	def print_details(self):

		print "num_faces = %d" % (self.num_faces)
		sleep(1)

		outline = ""
		for i in self.index:
			outline += "%d " % (i)
			
		print outline

	def write_to_file(self, fname):
		
		with open(fname, "w") as f:
			f.write("ffea vdw file\nnum_faces %d\nvdw params:\n" % (self.num_faces))
			for i in self.index:
				f.write("%d\n" % (i))

	def reset(self):

		self.index = []
		self.num_faces = 0
		self.valid = False
