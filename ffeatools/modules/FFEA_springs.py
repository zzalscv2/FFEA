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

from os import path
from time import sleep
import sys

class FFEA_springs:

	def __init__(self, fname = ""):
	
		self.reset()
		if fname == "":
			return

		try:
			self.load(fname)
		except:
			raise

	def load(self, fname):

		sys.stdout.write("Loading FFEA springs file...")
		
		# File format?
		base, ext = path.splitext(fname)
		try:
			if ext == ".spring" or ext == ".springs":

				self.load_springs(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".spring", ".springs"])

		except:
			raise

		sys.stdout.write("done!\n")


	def load_springs(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea springs file" and line != "walrus springs file":
			raise FFEAFormatError(lin=1, lstr="ffea springs file")

		try:
			self.num_springs = int(fin.readline().split()[1])

		except IndexError:
			raise FFEAFormatError(lin="2", lstr="num_springs %d")

		if fin.readline().strip() != "springs:":
			raise FFEAFormatError(lin="5", lstr="surface nodes:")

		# Read springs now
		try:
			for i in range(self.num_springs):
				sline = fin.readline().split()

				# Get a spring
				spring = FFEA_spring()
				spring.set_properties(sline[0], sline[1], sline[2:4], sline[4:6], sline[6:8])
				self.add_spring(spring)

		except (IndexError, ValueError):
			raise FFEAFormatError(lin=i+j+6, lstr="%f %f %d %d %d %d %d %d")
		except:
			raise

		fin.close()

	def add_spring(self, spring):

		self.spring.append(spring)
		self.num_springs += 1

	def get_num_springs(self):

		return self.num_springs

	def print_details(self):

		print "num_springs = %d" % (self.num_springs)
		sleep(1)

		print "\n\t k\t\tl\t\tblob\tconf\tnode"
		for s in self.spring:
			index = self.spring.index(s)
			outline = "Spring " + str(index) + " "
			outline += "%e\t%e\t%d %d\t%d %d\t%d %d" % (s.k, s.l, s.blob_index[0], s.blob_index[1], s.conformation_index[0], s.conformation_index[1], s.node_index[0], s.node_index[1])

			print outline

	def write_to_file(self, fname):
		
		fout = open(fname, "w")
		fout.write("ffea springs file\nnum_springs %d\nsprings:\n" % (self.num_springs))
		for s in self.spring:
			fout.write("%e %e %d %d %d %d %d %d\n" % (s.k, s.l, s.blob_index[0], s.blob_index[1], s.conformation_index[0], s.conformation_index[1], s.node_index[0], s.node_index[1]))
		fout.close()

	def reset(self):

		self.spring = []
		self.num_springs = 0

class FFEA_spring:

	def __init__(self):

		self.reset()

	def set_properties(self, k, l, bin, cin, nin):

		try:
			self.k = float(k)
			self.l = float(l)
			self.blob_index = [int(i) for i in bin]
			self.conformation_index = [int(i) for i in cin]
			self.node_index = [int(i) for i in nin]
		except:
			raise

	def print_details(self):
		print "\n\t k\t\tl\t\tblob\tconf\tnode"
		print "%e\t%e\t%d %d\t%d %d\t%d %d" % (self.k, self.l, self.blob_index[0], self.blob_index[1], self.conformation_index[0], self.conformation_index[1], self.node_index[0], self.node_index[1])

	def reset(self):
		
		self.k = 0
		self.l = 0
		self.blob_index = []
		self.conformation_index = []
		self.node_index = []
