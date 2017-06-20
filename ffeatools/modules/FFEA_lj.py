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

class FFEA_lj:

	def __init__(self, fname = ""):
	
		self.reset()

		# Return empty object if fname not initialised
		if fname == "":
			self.valid = True
			sys.stdout.write("done! Empty object initialised.\n")
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

		sys.stdout.write("Loading FFEA LJ file...")
	
		# File format?
		base, ext = os.path.splitext(fname)
		try:
			if ext == ".lj":
				self.load_lj(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".lj"])

		except:
			raise

		self.valid = True
		self.empty = False
		sys.stdout.write("done!\n")

	
	def load_lj(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found. Returning empty object...")
			self.reset()

		# Test format
		line = fin.readline().strip()
		if line != "ffea vdw forcefield params file" and line != "walrus vdw forcefield params file":
			print("\tExpected 'ffea vdw file' but found " + line)
			raise TypeError

		self.num_face_types = int(fin.readline().split()[1])


		# Read interaction matrix now
		for i in range(self.num_face_types):
			sline = fin.readline().split()
			self.interaction.append([])
			for s in sline:
				intline = s.strip()[1:-1].split(",")
				self.interaction[-1].append(FFEA_lj_pair())
				self.interaction[-1][-1].eps = float(intline[0])
				self.interaction[-1][-1].r = float(intline[0])

		fin.close()

	def default(self):

		# Default to steric (ish) params only
		self.interaction = [[FFEA_lj_pair() for i in range(self.num_face_types)] for j in range(self.num_face_types)]
		for i in range(self.num_face_types):
			for j in range(i, self.num_face_types):
				self.interaction[i][j].eps = 1e-15
				self.interaction[i][j].r = 1e-10

				# Symmetric
				self.interaction[j][i] = self.interaction[i][j]

	def set_interaction_pair(self, i, j, eps, r0):
		self.interaction[i][j].eps = eps
		self.interaction[i][j].r = r0

		# Symmetric
		self.interaction[j][i] = self.interaction[i][j]

	def write_to_file(self, fname):

		fout = open(fname, "w")
		fout.write("ffea vdw forcefield params file\nnum_vdw_face_types %d\n" % (self.num_face_types))
		for i in self.interaction:
			for j in i:
				fout.write("(%e, %e) " % (j.eps, j.r))
			fout.write("\n")
		fout.close()

	def reset(self):

		self.num_face_types = 7
		self.interaction = []
		self.default()

		self.valid = False
		self.empty = True

class FFEA_lj_pair:

	def __init__(self):

		self.reset()

	def reset(self):

		self.r = 0.0
		self.eps = 0.0
		
