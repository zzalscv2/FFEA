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

import os, sys
from FFEA_exceptions import *

class FFEA_skeleton:

	def __init__(self, fname = ""):
	
		self.reset()

		if fname == "":
			self.valid = True
			sys.stdout.write("Empty skeleton object initialised.\n")
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

		sys.stdout.write("Loading FFEA skeleton file...")

		# File format?
		base, ext = os.path.splitext(fname)
		try:
			if ext == ".skel":
				self.load_skel(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".skel"])

		except:
			raise

		self.valid = True
		self.empty = False
		sys.stdout.write("done!\n")

	def load_skel(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Test format
		line = fin.readline().strip()
		print(line)
		if line != "ffea skeleton file":
			print("\tExpected 'ffea skeleton file' but found " + line)
			raise TypeError

		# Object list sizes
		for i in range(2):
			sline = fin.readline().split()
			if("joint" in sline[0]):
				self.num_joints = int(sline[1])
			elif("bone" in sline[0]):
				self.num_bones = int(sline[1])
	
		# Read both sets of data
		start = fin.tell()
		while(True):
			line = fin.readline()
			if("joints:" in line):
				self.read_joints(fin)
				break
		fin.seek(start)
		while(True):
			line = fin.readline()
			if("bones:" in line):
				self.read_bones(fin)
				break

	def read_joints(self, fin):
		
		# Get a joints list
		self.joints = [0 for i in range(self.num_joints)]
		for i in range(self.num_joints):
			self.joints[i] = int(fin.readline())

	def read_bones(self, fin):
		
		# Get a joints list
		self.bones = [[0,0] for i in range(self.num_bones)]
		for i in range(self.num_bones):
			sline = fin.readline().split()
			for j in range(2):
				self.bones[i][j] = int(sline[j])

	def reset(self):

		self.num_joints = 0
		self.num_bones = 0
		self.bones = []
		self.joints = []
		self.valid = False
		self.empty = True
