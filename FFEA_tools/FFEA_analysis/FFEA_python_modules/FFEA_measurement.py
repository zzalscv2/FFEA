import sys
from os import path
import numpy as np
try:
  import matplotlib.pyplot as plt
except:
  plt = False

class FFEA_measurement:

	def __init__(self, fname):

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

			dfname = path.splitext(fname)[0] + ".fdm"
			if path.exists(fname):
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
		
		line = fin.readline().strip()
		while(line != "Measurements:"):
			if "num_blobs" in line:
				self.num_blobs = int(line.split("=")[1])

			line = fin.readline().strip()

		# Build a dictionary of possible variables and a map to those
		self.global_meas = {'Time': None, 'KineticEnergy': None, 'StrainEnergy': None, 'SpringEnergy': None, 'VdWEnergy': 'None', 'PreCompEnergy': None, 'Centroid.x': None, 'Centroid.y': None, 'Centroid.z': None, 'Centroid': None, 'RMSD': None}
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
		
		while(line.strip() != ""):
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
				self.blob_meas[i][key] = np.array(self.blob_meas[i][key])

			for j in range(i, self.num_blobs):
				for key in self.interblob_meas[i][j]:
					self.interblob_meas[i][j][key] = np.array(self.interblob_meas[i][j][key])
				self.interblob_meas[j][i] = self.interblob_meas[i][j]


	def reset(self):

		self.num_blobs = 0
		self.global_meas = None
		self.blob_meas = []
		self.interblob_meas = []
				

'''
class FFEA_measurement_blob:

	def __init__(self, fname = "", frame_rate=1, num_frames_to_read=1000000):
	
		self.reset()
		if fname == "":
			return

		try:
			self.load(fname, frame_rate=frame_rate, num_frames_to_read=num_frames_to_read)
		except:
			return
	
	def load(self, fname, frame_rate=1, num_frames_to_read=1000000):

		# Test file exists
		if not os.path.exists(fname):
			print("\tFile '" + fname + "' not found.")
			return

		# Read (dictionary and measurement keys should be made to line up in the future, or make so that we can have arbitrary extra measurements after the core ones isn't needed)
		try:
			with open(fname, "r") as fin:

				# Get to start
				while(fin.readline()[0]) != "#":
					pass

				# Then, assign stuff
				frame = -1
				line = fin.readline()
				while(line != ""):
					
					frame += 1

					# Should we end?
					if frame > num_frames_to_read:
						break

					# Should we skip?
					if frame % frame_rate != 0:
						line = fin.readline()
						continue

					sline = line.split()
					try:
						self.params['Step'].append(int(sline[0]))
						self.params['Kinetic'].append(float(sline[1]))
						self.params['Elastic'].append(float(sline[2]))
						self.params['Centre'].append([float(sline[3]), float(sline[4]), float(sline[5])])
						self.params['AMomentum'].append([float(sline[6]), float(sline[7]), float(sline[8])])
						self.params['RMSD'].append(float(sline[9]))
				
						line = fin.readline()
					except:
						line = fin.readline()
						continue
		
				# Make numpy for plotting and analysing and stuff
				for key in self.params:
					self.params[key] = np.array(self.params[key])
		except:
			print("\tUnable to load FFEA_measurement from " + fname + ". Returning empty object...")


	def plot(self, key, **opts):
		
		if self.params.get(key) == None:

			print("Error. '" + key + "' is not a recognised variable. Cannot be plotted.")
			return
		
		# Send to the appropriate function
		if key == "Elastic" or key == "Kinetic":
			self.plot_energy(opts, energy=key)
		
		else:
			print key + " is currently unsupported. Sorry. If you're a dev, program this in!"

		
	
	def plot_energy(self, opts, energy=None):
		
		print("Currently unavailable")
		return
		if plt == False: 
			print("plot_energy unavailable: matplotlib could not be imported")

		"""
		if energy == None:
			print "Error. 'energy' should be 'Elastic' or 'Kinetic'"
			return		
		
		# Default graphing options
		plotType = {"hist": False, "logx": False, "logy": False, "ravg": False}
		ylabel = energy + " Energy (J)"
		y = self.params[energy]
		print y
		savefig = [False, None]

		xlabel = "Step"
		x = self.params['Step']
		print x
		# Get opts (this is super lazy, and will make the later options take priority rather than raise an error. But! It makes them mutually exclusive)
		for o in opts:
			
			for t in plotType:
				if t == o:
					plotType[t] = True
				else:
					plotType[t] = False

		
		for o, val in opts.items():
			if o.lower() == "kt" or o.lower() == "kbt":
				ylabel = ylabel[0:-4] + "(kT)"
				y /= float(val)
			
			elif o.lower() == "dt" or o.lower() == "step":
				xlabel = "Time (ns)"
				x *= float(val)
				x /= 1e-9

			elif o == "save" or o == "savefig":
				savefig[0] = "True"
				savefig[1] = val

		# Do stuff
		if plotType['ravg'] == True:
			
			avg = []
			for i in range(len(y)):
				avg.append(np.mean(y[0:i]))

			plt.plot(x, avg)
			plt.title(energy + " Energy Running Time Average")

		elif plotType['logx'] == True:
			
			plt.semilogx(x, y)
			plt.title(energy + " Energy Trace (log x)")
			xlabel = "log " + xlabel[0:-4]

		elif plotType['logy'] == True:
			
			plt.semilogx(x, y)
			plt.title(energy + " Energy Trace (log y)")
			ylabel = "log " + ylabel[0:-4]

		elif plotType['hist'] == True:
			
			pass
		else:

			# Standard time trace
			plt.plot(x,y)
			plt.title(energy + " Energy Trace")

		ax = plt.gca()
		ax.set_ylim(bottom=0, top = max(y) * 1.5)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
	
		if savefig[0] == True:
			plt.savefig(savefig[1])

		plt.show()
		"""

	def reset(self):

		self.params = {'Step': [], 'Kinetic': [], 'Elastic': [], 'Centre': [], 'AMomentum': [], 'RMSD': []}

'''
