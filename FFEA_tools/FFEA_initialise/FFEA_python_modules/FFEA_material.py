import numpy as np
import sys

class FFEA_material:

	def __init__(self, fname):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. File " + fname  + " not found."
			return

		# Header
		if fin.readline().rstrip() != "ffea material params file":
			print "Error. Expected to read 'ffea material params file'. This may not be an ffea material params file"
			return

		# num_elements
		try:
			self.num_elements = int(fin.readline().split()[1])
			self.element = []

		except(ValueError):
			print "Error. Expected to read:"
			print "num_elements = %d"
			self.reset()
			fin.close()
			return

		# material parameters
		for i in range(self.num_elements):
			try:
				line = fin.readline()
				if line == [] or line == None or line == "":
					raise EOFError

				sline = line.split()
				self.element.append(FFEA_material_element(float(sline[0]), float(sline[1]), float(sline[2]), float(sline[3]), float(sline[4]), float(sline[5])))

			except(EOFError):
				print "Error. EOF may have been reached prematurely:\nnum_elements = " + str(self.num_elements) + "\nnum_elements read = " + str(i)
				self.reset()
				fin.close()
				return

			except(IndexError, ValueError):
				print "Error. Expected a material radius of the form '%f' for node " + str(i) + ", but found " + line
				self.reset()
				fin.close()
				return
		
		fin.close()

	def reset(self):
		self.element = []
		self.num_elements = 0

class FFEA_material_element:

	def __init__(self, d, sv, bv, sm, bm, de):

		self.density = d
		self.shear_viscosity = sv
		self.bulk_viscosity = bv
		self.shear_modulus = sm
		self.bulk_modulus = bm
		self.dieletric = de
