from os import path
from time import sleep
import numpy as np
from FFEA_topology import FFEA_topology
from FFEA_exceptions import *

class FFEA_material:

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

		print("Loading FFEA material file...")

		# Test file exists
		if not path.exists(fname):
			print fname + " not found."
			raise IOError
	
		# File format?
		base, ext = path.splitext(fname)
		if ext == ".mat":
			try:
				self.load_mat(fname)
				self.valid = True
			except:
				raise IOError

		else:
			raise FFEAIOError(fname=fname, fext=[".mat"])

	def load_mat(self, fname):

		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			raise

		# Test format
		line = fin.readline().strip()
		if line != "ffea material params file" and line != "walrus material params file":
			raise FFEAFormatError(lin=1, lstr="ffea material params file")

		try:
			num_elements = int(fin.readline().split()[1])
		except:
			raise FFEAFormatError(lin="2", lstr="num_elements %d\n")

		# Read elements now
		try:
			for i in range(num_elements):
				sline = fin.readline().split()
	
				# Get an element (we want a matrix of values for slicing)
				el = []

				for s in sline:
					el.append(float(s))
				self.add_element(el)

		except (IndexError, ValueError):
			raise FFEAFormatError(lin=i+3, lstr="%f %f %f %f %f %f")
		except:
			raise

		fin.close()

		# Numpy the params up, as they are basically a matrix
		self.element = np.array(self.element)

	def build(self, num_elements, **params):
		
		# Get a number of elements
		if num_elements == 0:
			self.reset()
			return
		
		# Check params are valid
		plist = []
		order = ["d", "sv", "bv", "sm", "bm", "di"]
		for key in order:
			try:
				plist.append(float(params[key]))
			except(KeyError):
				print "Error. Required build parameters not present. We need:"
				for key2 in order:
					print "\t" + key2

				raise
		
		# Get params
		# Now build in correct way
		for i in range(num_elements):
			self.add_element(plist)
	
		# Numpy the params up, as they are basically a matrix
		self.element = np.array(self.element)

	def write_to_file(self, fname):
	
		fout = open(fname, "w")
		fout.write("ffea material params file\nnum_elements %d\n" % (self.num_elements))
		for el in self.element:
			fout.write("%6.3f %6.3f %6.3f %10.1f %10.1f %6.3f" % (el[0], el[1], el[2], el[3], el[4], el[5]))
			fout.write("\n")
		fout.close()
		
	def add_element(self, el):

		self.element.append(el)
		self.num_elements += 1

	def set_params(self, index, d, sv, bv, sm, bm, di):
	
		try:
			self.element[index] = [float(d), float(sv), float(bv), float(sm), float(bm), float(di)]
		except(IndexError):
			print "Element " + str(index) + " does not yet exist."
			 
	def get_num_elements(self):

		return len(self.element)

	def print_details(self):

		print "num_elements = %d" % (self.num_elements)

		print "\t\tDensity,Shear Viscosity,Bulk Viscosity,Shear Modulus,Bulk Modulus,Dielectric Constant\n"
		sleep(1)

		index = -1
		for e in self.element:
			index += 1
			outline = "Element " + str(index) + "\t"
			for param in e:
				outline += str(param) + ", "

			print outline

	def reset(self):

		self.element = []
		self.num_elements = 0
		self.valid = False 
