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
import numpy as np

class FFEA_kinetic_states:

	def __init__(self, fname = "", frame = 0):
	
		self.reset()

		try:
			self.load(fname)
		except:
			return

	def load(self, fname, findex = 0):

		print("Loading FFEA kinetic states file...")

		# Test file exists
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".states":
			try:
				self.load_states(fname)
			except:
				print("\tUnable to load FFEA_kinetic_states from " + fname + ". Returning empty object...")

		else:
			print("\tUnrecognised file extension '" + ext + "'.")

	def load_states(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea kinetic states file" and line != "walrus kinetic states file":
			print("\tExpected 'ffea kinetic states file' but found " + line)
			raise TypeError

		num_states = int(fin.readline().split()[1])

		fin.readline()

		# Read states now
		while(True):
			try:
				conformation_index = int(fin.readline().split()[1])
				sites = fin.readline().split()[1:]
				
			except(IndexError):
				break

			# Get a state
			state = FFEA_kinetic_state()
			state.set_conformation(conformation_index)
			state.set_sites(sites)
			self.state.append(state)

		fin.close()
		self.num_states = len(self.state)
		
	def print_details(self):

		print("num_states = %d" % (self.num_states))
		sleep(1)


		print("states:\n")

		index = -1
		for s in self.state:

			index += 1
			outline = ""
			outline += "State " + str(index) + ":\n\tConformation index = " + str(s.conformation_index)
			if s.bound:
				outline += "\n\tBound from site type " + str(s.site[0]) + " to type " + str(s.site[1]) + "\n"

			print(outline)
	
	def reset(self):

		self.num_states = 0
		self.state = []

class FFEA_kinetic_state:

	def __init__(self):
		
		self.reset()

	def set_conformation(self, index):
		self.conformation_index = index

	def set_sites(self, alist):
		
		try:
			self.site = [int(alist[0]), int(alist[1])]
		except(IndexError):
			return

		self.bound = True

	def reset(self):

		self.bound = False
		self.site = [-1,-1]
		self.conformation_index = 0
		
