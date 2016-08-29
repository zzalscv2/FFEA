import sys
from os import path
import numpy as np
try:
  import matplotlib.pyplot as plt
except:
  plt = False

class FFEA_measurement:

	def __init__(self, fname = ""):

		self.reset()

		# Test what type of file it is
		if not path.exists(fname):
			print("\tFile '" + fname + "' not found.")
			return
		
		fin = open(fname, "r")
		line = fin.readline().strip()
		fin.close()
		try:
			if line == "FFEA Global Measurement File":
				self.load_global(fname)
			else:		
				print("\tPlease supply us with the global measurement file, not the '-d' .fdm file")				
				return

			# Get num frames for quick access
			self.num_frames = len(self.global_meas["Time"])

			dfname = path.splitext(fname)[0] + ".fdm"
			if path.exists(dfname):
				self.load_detailed(dfname)
			
		except:
			self.reset()
			return
	
	def load_global(self, fname):

		print("Loading FFEA Global Measurement file...")
	
		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise
		
		# Details first
		line = fin.readline().strip()
		while(line != "Simulation Details:"):	
			line = fin.readline().strip()

		# Sim time
		line = fin.readline().split()
		self.date = [int(i) for i in line[3].split("/")]
		self.time = [int(i) for i in line[5].split(":")]

		# Script fname
		line = fin.readline().strip()
		line = line.split("=")
		if len(line) != 1:
			self.script_fname = line[1].strip()
		else:
			self.script_fname = ""

		# Sim type
		line = fin.readline().strip()
		line = line.split("=")
		if len(line) != 1:
			self.simtype = line[1].strip()
		else:
			self.simtype = "Full"

		line = fin.readline().strip()
		while(line != "Measurements:"):
			if "num_blobs" in line:
				self.num_blobs = int(line.split("=")[1])

			line = fin.readline().strip()

		# Build a dictionary of possible variables and a map to those
		self.global_meas = {'Time': None, 'KineticEnergy': None, 'StrainEnergy': None, 'SpringEnergy': None, 'VdWEnergy': None, 'PreCompEnergy': None, 'Centroid.x': None, 'Centroid.y': None, 'Centroid.z': None, 'Centroid': None, 'RMSD': None}
		sline = fin.readline().strip().split()
		measmap = ["" for i in range(len(sline))]		
		i = -1
		for title in sline:
			i += 1
			measmap[i] = title.strip()
			self.global_meas[measmap[i]] = []

		# Read measurements
		line = fin.readline()
		while(line != ""):
			if line.strip() == "#==RESTART==":
				line = fin.readline()
				continue

			sline = line.split()
			for i in range(len(sline)):
				self.global_meas[measmap[i]].append(float(sline[i]))
			line = fin.readline()
			
		# Move centroid into more useful format
		self.global_meas["Centroid"] = []
		for i in range(len(self.global_meas['Time'])):
			self.global_meas["Centroid"].append([self.global_meas["Centroid.x"][i], self.global_meas["Centroid.y"][i], self.global_meas["Centroid.z"][i]])

		del self.global_meas["Centroid.x"]
		del self.global_meas["Centroid.y"]
		del self.global_meas["Centroid.z"]
	
		for key in self.global_meas:
			if self.global_meas[key] != None:
				self.global_meas[key] = np.array(self.global_meas[key])

	def load_detailed(self, fname):

		print("Loading FFEA Detailed Measurement file...")
	
		# Open file
		try:
			fin = open(fname, "r")
		except(IOError):
			print("\tFile '" + fname + "' not found.")
			self.reset()
			raise

		line = fin.readline().strip()
		while(line != "Measurements:"):
			line = fin.readline().strip()

		# Get column title line and build maps to the required stuff
		line = fin.readline()
		sline = line.split("|")[1:]
		localsline = sline[:self.num_blobs]
		globalsline = sline[self.num_blobs:]
		
		# Build dictionaries ans maps to variables
		
		# Local to blobs first
		# Initialise arrays (measmap and indexmap are of unknown length initially)
		self.blob_meas = [{'KineticEnergy': None, 'StrainEnergy': None, 'Centroid.x': None, 'Centroid.y': None, 'Centroid.z': None, 'Centroid': None, 'RMSD': None} for i in range(self.num_blobs)]
		indexmap = []
		measmap = []

		for i in range(self.num_blobs):
			ssline = localsline[i].split()[1:]
			for title in ssline:
				indexmap.append(i)
				measmap.append(title.strip())
				self.blob_meas[i][title.strip()] = []

		# Now the global one
		# Initialise arrays
		self.interblob_meas = [[{"VdWEnergy": None, "SpringEnergy": None, "PreCompEnergy": None} for i in range(self.num_blobs)] for j in range(self.num_blobs)]
		iindexmap = []
		imeasmap = []

		for i in range(len(globalsline)):

			# Get the pair of indices
			ssline = globalsline[i].split()
			indexpair = [int(j) for j in ssline[0][1:].split("B")]
			ssline = ssline[1:]
			
			for title in ssline:
				iindexmap.append(indexpair)
				imeasmap.append(title.strip())
				self.interblob_meas[indexpair[0]][indexpair[1]][title.strip()] = []

		# Now, read measurements and fill the relevent arrays
		line = fin.readline()
		
		while(line != ""):
			if line.strip() == "#==RESTART==":
				line = fin.readline()
				continue

			sline = line.split()[1:]

			localsline = sline[:len(indexmap)]
			globalsline = sline[len(indexmap):]

			# Local to blobs first!
			for i in range(len(localsline)):
				self.blob_meas[indexmap[i]][measmap[i]].append(float(localsline[i]))

			# Now global
			for i in range(len(globalsline)):
				self.interblob_meas[iindexmap[i][0]][iindexmap[i][1]][imeasmap[i]].append(float(globalsline[i]))

			line = fin.readline()
		

		# Move centroid into more useful format, make interblob array symmetric and turn stuff to numpy
		for i in range(self.num_blobs):
			self.blob_meas[i]["Centroid"] = []

			for j in range(len(self.global_meas['Time'])):
				self.blob_meas[i]["Centroid"].append([self.blob_meas[i]["Centroid.x"][j], self.blob_meas[i]["Centroid.y"][j], self.blob_meas[i]["Centroid.z"][j]])

			del self.blob_meas[i]["Centroid.x"]
			del self.blob_meas[i]["Centroid.y"]
			del self.blob_meas[i]["Centroid.z"]

			for key in self.blob_meas[i]:
				if self.blob_meas[i][key] != None:
					self.blob_meas[i][key] = np.array(self.blob_meas[i][key])

			for j in range(i, self.num_blobs):
				for key in self.interblob_meas[i][j]:
					if self.interblob_meas[i][j][key] != None:
						self.interblob_meas[i][j][key] = np.array(self.interblob_meas[i][j][key])
				self.interblob_meas[j][i] = self.interblob_meas[i][j]


	def add_empty_blob(self):
		
		if self.global_meas == None:
			self.global_meas = self.global_meas = {'Time': None, 'KineticEnergy': None, 'StrainEnergy': None, 'SpringEnergy': None, 'VdWEnergy': None, 'PreCompEnergy': None, 'Centroid.x': None, 'Centroid.y': None, 'Centroid.z': None, 'Centroid': None, 'RMSD': None}
		self.blob_meas.append({'KineticEnergy': None, 'StrainEnergy': None, 'Centroid.x': None, 'Centroid.y': None, 'Centroid.z': None, 'Centroid': None, 'RMSD': None})
		for i in range(self.num_blobs):
			self.interblob_meas[i].append([{"VdWEnergy": None, "SpringEnergy": None, "PreCompEnergy": None}])
		self.interblob_meas.append([{"VdWEnergy": None, "SpringEnergy": None, "PreCompEnergy": None} for i in range(self.num_blobs + 1)])
		self.num_blobs += 1

	def write_to_file(self, fname, script = None):

		# Get filenames first
		globalfname = fname
		blobfname = path.splitext(fname)[0] + ".fdm"

		print "Writing Measurements to file:\n\tGlobal Will be written to %s" % (globalfname)
		if self.blob_meas != []:
			print "\tDetailed will be written to %s\n" % (blobfname)


		#
		# Global first
		#
		
		fout = open(globalfname, "w")
		fout.write("FFEA Global Measurement File\n\nSimulation Details:\n")
		fout.write("Simulation Began on %d/%d/%d at %d:%d:%d\n" % (self.date[0],self.date[1], self.date[2], self.time[0], self.time[1], self.time[2]))
		fout.write("Script Filename = %s" % (self.script_fname))
		fout.write("Simulation Type = %s\n\n" % (self.simtype))

		# Params, maybe		
		if script != None:
			script.params.write_to_file(fout, self.script_fname)

		# Measurements
		keys_to_write = ["Time"]
		if self.global_meas["KineticEnergy"] != None:
			keys_to_write.append("KineticEnergy")
		keys_to_write.append("StrainEnergy")
		if self.global_meas["SpringEnergy"] != None:
			keys_to_write.append("SpringEnergy")
		if self.global_meas["VdWEnergy"] != None:
			keys_to_write.append("VdWEnergy")
		if self.global_meas["PreCompEnergy"] != None:
			keys_to_write.append("PreCompEnergy")
		keys_to_write.append("Centroid")
		keys_to_write.append("RMSD")
		fout.write("\nMeasurements:\n")
		for key in keys_to_write:
			if key == "Centroid":
				fout.write("%-14s%-14s%-14s" % ("Centroid.x", "Centroid.y", "Centroid.z"))
			else:
				fout.write("%-14s" % (key))
		fout.write("\n")
		for i in range(self.num_frames):
			for key in keys_to_write:
				if key == "Centroid":
					fout.write("%-14.6e%-14.6e%-14.6e" % (self.global_meas[key][i][0],self.global_meas[key][i][1], self.global_meas[key][i][2]))
				else:
					fout.write("%-14.6e" % (self.global_meas[key][i]))
			fout.write("\n")
		fout.close()
		
		#
		# Now detailed bit
		#

		if self.blob_meas != []:
			fout = open(blobfname, "w")
			fout.write("FFEA Detailed Measurement File\n\nMeasurements:\n")

			# We need a keys_to_write for each blob
			keys_to_write = [[] for i in range(self.num_blobs)]
			pair_keys_to_write = [[[] for j in range(i, self.num_blobs)] for i in range(self.num_blobs)]
			do_interblob = [[False for j in range(i, self.num_blobs)] for i in range(self.num_blobs)]
			for i in range(self.num_blobs):
				if self.blob_meas[i]["KineticEnergy"] != None:
					keys_to_write[i].append("KineticEnergy")
				keys_to_write[i].append("StrainEnergy")
				keys_to_write[i].append("Centroid")
				keys_to_write[i].append("RMSD")
				for j in range(i, self.num_blobs):
					if self.interblob_meas[i][j]["VdWEnergy"] != None:
						pair_keys_to_write[i][j].append("VdWEnergy")
						do_interblob[i][j] = True
					if self.interblob_meas[i][j]["SpringEnergy"] != None:
						pair_keys_to_write[i][j].append("SpringEnergy")
						do_interblob[i][j] = True
					if self.interblob_meas[i][j]["PreCompEnergy"] != None:
						pair_keys_to_write[i][j].append("PreCompEnergy")
						do_interblob[i][j] = True

		

			# Time first
			fout.write("%-14s" % ("Time"))
		
			# Then local to blobs
			for i in range(self.num_blobs):
				fout.write("| B%d " % (i))
				for key in keys_to_write[i]:
					if key == "Centroid":
						fout.write("%-14s%-14s%-14s" % ("Centroid.x", "Centroid.y", "Centroid.z"))
					else:
						fout.write("%-14s" % (key))

			# Then interblob stuff
			if do_interblob:
				for i in range(self.num_blobs):
					for j in range(i, self.num_blobs):
						if do_interblob[i][j]:
							fout.write("| B%dB%d " % (i, j))
							for key in pair_keys_to_write[i][j]:
								fout.write("%-14s" % (key))
			fout.write("\n")

			for i in range(self.num_frames):
				fout.write("%-14.6e" % (self.global_meas["Time"][i]))
				for j in range(self.num_blobs):
					fout.write("     ")
					for key in keys_to_write[j]:
						if key == "Centroid":
							fout.write("%-14.6e%-14.6e%-14.6e" % (self.blob_meas[j][key][i][0],self.blob_meas[j][key][i][1], self.blob_meas[j][key][i][2]))
						else:
							fout.write("%-14.6e" % (self.blob_meas[j][key][i]))

				for j in range(self.num_blobs):
					for k in range(j, self.num_blobs):
						if do_interblob[j][k]:
							fout.write("      ")
							for key in pair_keys_to_write[j][k]:
								fout.write("%-14.6e" % (self.interblob_meas[j][k][key][i]))
				fout.write("\n")

			fout.close()

	def reset(self):

		self.date = None
		self.time = None
		self.script_fname = ""
		self.num_blobs = 0
		self.num_frames = 0
		self.global_meas = None
		self.blob_meas = []
		self.interblob_meas = []
