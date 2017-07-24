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
import copy
import numpy as np
from FFEA_exceptions import *

class FFEA_ctforce_linear:

	def __init__(self):

		self.reset()

	def set_params(self, fmag, fvec, bin, cin, nin):

		self.force_mag = fmag
		self.force_vector = np.array(fvec)
		self.bindex = bin
		self.cindex = cin
		self.nindex = nin

	def reset(self):

		self.force_mag = 0.0
		self.force_vector = np.array([0.0,0.0,0.0])
		self.bindex = 0
		self.cindex = 0
		self.nindex = 0

class FFEA_ctforce_rotational:

	def __init__(self):

		self.reset()

	def set_params(self, fmag, apoint, axis, bin, cin, nin, kw):

		self.force_mag = fmag
		self.axis_point = np.array(apoint)
		self.axis = np.array(axis)
		self.bindex = bin
		self.cindex = cin
		self.nindex = nin
		self.keyword = kw

	def reset(self):

		self.force_mag = 0.0
		self.axis_point = np.array([0.0,0.0,0.0])
		self.axis = np.array([0.0,0.0,0.0])
		self.bindex = 0
		self.cindex = 0
		self.nindex = 0
		self.keyword = None

class FFEA_ctforce_linear_surface:

	def __init__(self):

		self.reset()

	def set_params(self, fmag, fvec, bin, cin, flist):

		self.force_mag = fmag
		self.force_vector = np.array(fvec)
		self.bindex = bin
		self.cindex = cin
		self.findex = flist
		self.num_faces = len(self.findex)

	def reset(self):

		self.force_mag = 0.0
		self.force_vector = np.array([0.0,0.0,0.0])
		self.bindex = 0
		self.cindex = 0
		self.num_faces = 0
		self.findex = []

class FFEA_ctforces:

	def __init__(self, fname = ""):
	
		self.reset()

		# Empty fname give an empty object
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

		sys.stdout.write("Loading FFEA CTForces file...")
	
		# File format?
		base, ext = os.path.splitext(fname)
		try:
			if ext == ".ctforces":

				self.load_ctforces(fname)
			else:
				raise FFEAIOError(fname=fname, fext=[".ctforces"])

		except:
			raise

		self.valid = True
		self.empty = False
		sys.stdout.write("done!\n")

	def load_ctforces(self, fname):

		
		# Open file
		try:
			fin = open(fname, "r")

		except(IOError):
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea ctforces file":
			raise FFEAFormatError(lin=1, lstr="ffea ctforces file")

		# Get all forces
		num_ctforces = copy.deepcopy(self.num_ctforces)
		for i in range(4):
			sline = fin.readline().split()
			try:
				if sline[0].strip() == "num_ctforces":
					num_ctforces["all"] = int(sline[1])

				elif sline[0].strip() == "num_linear_forces":
					num_ctforces["lin"] = int(sline[1])

				elif sline[0].strip() == "num_rot_forces":
					num_ctforces["rot"] = int(sline[1])

				elif sline[0].strip() == "num_linear_surface_forces":
					num_ctforces["lin_surf"] = int(sline[1])

				else:
					raise IndexError

			except:
				raise FFEAFormatError(lin="2-5", lstr="num_ctforces %d\nnum_linear_forces %d\nnum_rot_forces %d\nnum_linear_surface_forces %d")

		# Test for consistency
		asum = 0
		for ftype in self.force_types[1:]:
			asum += num_ctforces[ftype]

		if asum != num_ctforces["all"]:
			print("Error. num_ctforces %d is inconsistent with forces specified in header.")			
			raise FFEAFormatError(lin="2", lstr="num_ctforces %d")

		# Get remaining lines
		force_line = fin.readlines()

		# Use to read forces
		if num_ctforces["lin"] != 0:
			sin = 0
			stline = force_line[sin].strip()
			while stline != "linear forces:":
				sin += 1
				try:
					stline = force_line[sin].strip()
				except(IndexError):
					raise FFEAFormatError(lin="", lstr="linear forces:")					

			for i in range(sin + 1, sin + 1 + num_ctforces["lin"]):
				sline = force_line[i].split()
				try:
					self.add_linear_force(sline[0], [sline[1], sline[2], sline[3]], sline[4], sline[5], sline[6])
				except:
					raise FFEAFormatError(lin="", lstr="force (mag) force (vec) blob conf node\n%f %f %f %f %d %d %d")
			
		if num_ctforces["rot"] != 0:
			sin = 0
			stline = force_line[sin].strip()
			while stline != "rotational forces:":
				sin += 1
				try:
					stline = force_line[sin].strip()
				except(IndexError):
					raise FFEAFormatError(lin="", lstr="rotational forces:")

			for i in range(sin + 1, sin + 1 + num_ctforces["rot"]):
				sline = force_line[i].split()
				
				# Check type keyword				
				if sline[1].strip()[0] == "n":
					print("Error. 'n' type rotational force not yet supported.")			
					raise FFEAFormatError(lin="", lstr="force (mag) keyword axis (point vec) axis (direction vec) blob conf node\n%f pt %f %f %f %d %d %d %d %d %d")

				try:
					self.add_rotational_force(sline[0], sline[1], [sline[2], sline[3],sline[4]], [sline[5], sline[6], sline[7]], sline[8], sline[9], sline[10])
				except:
					raise FFEAFormatError(lin="", lstr="force (mag) keyword axis (point vec) axis (direction vec) blob conf node\n%f pt %f %f %f %f %f %f %d %d %d")

		if num_ctforces["lin_surf"] != 0:
			sin = 0
			stline = force_line[sin].strip()
			while stline != "linear surface forces:":
				sin += 1
				try:
					stline = force_line[sin].strip()
				except(IndexError):
					raise FFEAFormatError(lin="", lstr="linear surface forces:")

			for i in range(sin + 1, sin + 1 + num_ctforces["lin_surf"]):
				sline = force_line[i].split()

				try:
					self.add_linear_surface_force(sline[0], [sline[1], sline[2], sline[3]], sline[4], sline[5], sline[6:])
				except:
					raise FFEAFormatError(lin="", lstr="force (mag) force (vec)  blob conf face_list\n%f %f %f %f %d %d   %d %d %d %d...")

	def add_linear_force(self, fmag, fvec, bin, cin, nin):

		# Let's clean up the forces at least
		fvec = np.array([float(f) for f in fvec])
		fvec *= 1.0 / np.linalg.norm(fvec)

		# Convert inputs
		try:
			fmag = float(fmag)
			bin = int(bin)
			cin = int(cin)
			nin = nin.strip()
			if nin != 'all':
				nin = int(nin)

		except(TypeError):
			raise

		# Test force params
		if fmag < 0:
			raise IndexError
		
		if bin < 0 or cin < 0 or nin < 0:
			raise IndexError

		# Make new linear force object and add to list
		l = FFEA_ctforce_linear()
		l.set_params(fmag, fvec, bin, cin, nin)
		self.ctforces["lin"].append(l)
		self.num_ctforces["lin"] += 1
		self.num_ctforces["all"] += 1

	def add_rotational_force(self, fmag, kw, apoint, axis, bin, cin, nin):

		# Let's clean up the forces at least
		apoint = np.array([float(f) for f in apoint])

		axis = np.array([float(a) for a in axis])
		axis *= 1.0 / np.linalg.norm(axis)

		# Convert inputs
		try:
			fmag = float(fmag)
			bin = int(bin)
			cin = int(cin)
			nin = nin.strip()
			if nin != 'all':
				nin = int(nin)

		except(TypeError):
			raise

		# Test force params
		if fmag < 0:
			raise IndexError
		
		if bin < 0 or cin < 0 or nin < 0:
			raise IndexError

		# Make new rotational force object and add to list
		r = FFEA_ctforce_rotational()
		r.set_params(fmag, apoint, axis, bin, cin, nin, kw)
		self.ctforces["rot"].append(r)
		self.num_ctforces["rot"] += 1
		self.num_ctforces["all"] += 1

	def add_linear_surface_force(self, fmag, fvec, bin, cin, flist):

		# Let's clean up the forces at least
		fvec = np.array([float(f) for f in fvec])
		fvec *= 1.0 / np.linalg.norm(fvec)

		# Convert inputs
		try:
			fmag = float(fmag)
			bin = int(bin)
			cin = int(cin)
			if 'all' not in flist:
				flist = [int(f) for f in flist]

		except(TypeError):
			raise

		# Test force params
		if fmag < 0:
			raise IndexError
		
		if bin < 0 or cin < 0:
			raise IndexError

		for f in flist:
			if f < 0:
				raise IndexError

		# Make new rotational force object and add to list
		ls = FFEA_ctforce_linear_surface()
		ls.set_params(fmag, fvec, bin, cin, flist)
		self.ctforces["lin_surf"].append(ls)
		self.num_ctforces["lin_surf"] += 1
		self.num_ctforces["all"] += 1

	def write_to_file(self, fname):

		print("Writing to " + fname + "...")
		fout = open(fname, "w")
		fout.write("ffea ctforces file\n")
		fout.write("num_ctforces %d\n" % (self.num_ctforces['all']))
		fout.write("num_linear_forces %d\n" % (self.num_ctforces['lin']))
		fout.write("num_rot_forces %d\n" % (self.num_ctforces['rot']))
		fout.write("num_linear_surface_forces %d\n" % (self.num_ctforces['lin_surf']))

		if(self.num_ctforces['lin'] != 0):
			fout.write("linear forces:\n")
			for f in self.ctforces['lin']:
				fout.write("%5.3e %3.2f %3.2f %3.2f %d %d %s\n" % (f.force_mag, f.force_vector[0], f.force_vector[1], f.force_vector[2], f.bindex, f.cindex, str(f.nindex)))
			fout.write("\n")

		if(self.num_ctforces['rot'] != 0):
			fout.write("rotational forces:\n")
			for f in self.ctforces['rot']:
				fout.write("%5.3e %s %5.3e %5.3e %5.3e %3.2f %3.2f %3.2f %d %d %s\n" % (f.force_mag, f.keyword, f.axis_point[0], f.axis_point[1], f.axis_point[2], f.axis[0], f.axis[1], f.axis[2], f.bindex, f.cindex, str(f.nindex)))
			fout.write("\n")

		if(self.num_ctforces['lin_surf'] != 0):
			fout.write("linear surface forces:\n")
			for f in self.ctforces['lin_surf']:
				fout.write("%5.3e %3.2f %3.2f %3.2f %d %d" % (f.force_mag, f.force_vector[0], f.force_vector[1], f.force_vector[2], f.bindex, f.cindex))
				if 'all' in f.findex:
					fout.write(" %s\n" % ('all'))
				else:
					for findex in f.findex:
						fout.write(" %d" % (findex))
					fout.write("\n")
			fout.write("\n")
		
		fout.close()
		print("done!")

	def reset(self):

		self.valid = False
		self.empty = True

		self.force_types = ["all", "lin", "rot", "lin_surf"]
		self.num_ctforces = {}
		for ftype in self.force_types:
			self.num_ctforces[ftype] = 0

		self.ctforces = {}
		for ftype in self.force_types:
			self.ctforces[ftype] = []
