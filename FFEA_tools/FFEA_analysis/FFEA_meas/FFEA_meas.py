import sys, os
from math import *
import numpy as np

class FFEA_meas:
		
	def __init__(self, meas_basename, num_blobs, num_frames_to_read):
		
		self.num_blobs = num_blobs;
		self.blob = []
		self.world = None
		self.num_frames = 0;
		for i in range(self.num_blobs):
			self.read_blob_meas_from_file(meas_basename, i, num_frames_to_read)

		self.read_world_meas_from_file(meas_basename, self.num_blobs, num_frames_to_read)

	def read_blob_meas_from_file(self, basename, index, num_frames):

		# Open file
		meas_fname = basename + "_blob" + str(index) + ".out"
		print "Reading FFEA measurement file " + meas_fname
		meas = open(meas_fname, "r")

		# Get initial crap
		#if(meas.readline().strip() != "FFEA_blob_measurement_file"):
			#sys.exit("Error. Expected 'FFEA_blob_measurement_file' line. This may not be an FFEA measurement file. Bye!")

		meas.readline()

		# Begin importing data
		print("\n\tReading measurement data")
		check = 0
		if num_frames <= 0:
			num_frames = 10000
		completed = 0

		ke = []
		pe = []
		comx = []
		comy = []
		comz = []
		angmomx = []
		angmomy = []
		angmomz = []
		rmsd = []
		vdw_area = []
		vdw_force = []
		vdw_energy = []

		for i in range(num_frames):
			
			# Exit if finished
			if completed == 1:
				break
			
			# Get a line
			line = meas.readline()

			# Check for eof
			if line == "":
				print "Specified " + str(num_frames) + " to read, but eof reached at " + str(i - 1) + "."
				print "Continuing...\n"	
				completed = 1
				check -= 1				
				break

			# Split line and assign stuff
			sline = line.split()
			if(float(sline[1]) < 0.0):
				ke.append(0.0)
			else:
				ke.append(float(sline[1]))

			if(float(sline[2]) < 0.0):
				pe.append(0.0)
			else:
				pe.append(float(sline[2]))

			comx.append(float(sline[3]))
			comy.append(float(sline[4]))
			comz.append(float(sline[5]))
			angmomx.append(float(sline[6]))
			angmomy.append(float(sline[7]))
			angmomz.append(float(sline[8]))
			rmsd.append(float(sline[9]))
			vdw_area.append(float(sline[10]))
			vdw_force.append(float(sline[11]))
			vdw_energy.append(float(sline[12]))

			check += 1
			if num_frames >= 100 and check % (num_frames / 50) == 0:
				print "\t\tRead " + str(check) + " frames"

		# Add a blob
		print "Done. Read " + str(check) + " frames in total."
		num_frames = check
		self.blob.append(FFEA_meas_blob(num_frames, np.array(ke), np.array(pe), np.array([comx, comy, comz]), np.array([angmomx, angmomy, angmomz]), np.array(rmsd), np.array(vdw_area), np.array(vdw_force), np.array(vdw_energy)))

	def read_world_meas_from_file(self, basename, num_blobs, num_frames):

		# Open file
		meas_fname = basename + "_world.out"
		print "Reading FFEA measurement file " + meas_fname
		meas = open(meas_fname, "r")

		# Get initial crap
		#if(meas.readline().strip() != "FFEA_blob_measurement_file"):
			#sys.exit("Error. Expected 'FFEA_blob_measurement_file' line. This may not be an FFEA measurement file. Bye!")

		meas.readline()

		# Begin importing data
		print("\n\tReading measurement data")
		check = 0
		if num_frames <= 0:
			num_frames = 10000
		completed = 0

		vdw_area = [[np.array([0.0 for j in range(num_frames)]) for i in range(num_blobs)] for k in range(num_blobs)]
		vdw_force = [[np.array([0.0 for j in range(num_frames)]) for i in range(num_blobs)] for k in range(num_blobs)]
		vdw_energy = [[np.array([0.0 for j in range(num_frames)]) for i in range(num_blobs)] for k in range(num_blobs)]
		for i in range(num_frames):
			
			# Exit if finished
			if completed == 1:
				break
			
			# Get a line
			line = meas.readline()

			# Check for eof
			if line == "":
				print "Specified " + str(num_frames) + " to read, but eof reached at " + str(i - 1) + "."
				print "Continuing...\n"	
				completed = 1
				check -= 1				
				break

			# Split line and assign stuff
			sline = line.split()
			sline_index = 1
			for j in range(num_blobs):
				vdw_area[j][j][i] = 0.0
				vdw_force[j][j][i] = 0.0
				vdw_energy[j][j][i] = 0.0
				for k in range(j + 1, num_blobs, 1):
					vdw_area[j][k][i] = float(sline[sline_index])
					vdw_force[j][k][i] = float(sline[sline_index + 1])
					vdw_energy[j][k][i] = float(sline[sline_index + 2])

					vdw_area[k][j][i] = vdw_area[j][k][i]
					vdw_force[k][j][i] = vdw_force[j][k][i]
					vdw_energy[k][j][i] = vdw_energy[j][k][i]

					sline_index += 3
		
			check += 1
			if num_frames >= 100 and check % (num_frames / 50) == 0:
				print "\t\tRead " + str(check) + " frames"
		
		# Add a world
		print "Done. Read " + str(check) + " frames in total."
		num_frames = check
		self.world = FFEA_meas_world(num_frames, num_blobs, vdw_area, vdw_force, vdw_energy)

	def calc_avg_energies(self):

		for blob in self.blob:
			blob.calc_avg_energies()

class FFEA_meas_blob:

	def __init__(self, num_frames, ke, pe, com, angmom, rmsd, vdw_area, vdw_force, vdw_energy):
		
		self.num_frames = num_frames
		self.ke = ke
		self.ke_avg = 0.0
		self.pe_avg = 0.0
		self.pe = pe
		self.com = com
		self.angmom = angmom
		self.rmsd = rmsd
		self.vdw_area = vdw_area
		self.vdw_force = vdw_force
		self.vdw_energy = vdw_energy

	def calc_avg_energies(self):
		
		self.ke_avg = np.mean(self.ke)
		self.pe_avg = np.mean(self.pe)

	def get_avg_energies(self):
		return self.ke_avg, self.pe_avg
		
class FFEA_meas_world:

	def __init__(self, num_frames, num_blobs, vdw_area, vdw_force, vdw_energy):

		self.num_frames = num_frames
		self.num_blobs = num_blobs
		self.vdw_area = vdw_area
		self.vdw_force = vdw_force
		self.vdw_energy = vdw_energy
