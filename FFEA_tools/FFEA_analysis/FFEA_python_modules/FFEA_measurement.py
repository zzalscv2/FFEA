import sys, os
import numpy as np
import matplotlib.pyplot as plt

class FFEA_measurement:

	def __init__(self, fname, num_blobs, num_frames_to_read = float("inf")):

		# Initialise stuff
		self.reset()

		self.num_blobs = num_blobs

		# Split name into world and blob names
		print fname
		basename, ext = os.path.splitext(os.path.abspath(fname))
		world_fname = basename + "_world" + ext
		blob_fname = [basename + "_blob" + str(i) + ext for i in range(self.num_blobs)]
		
		# Start reading the files
		for fname in blob_fname:
			ablob = FFEA_blob_measurement(fname, num_frames_to_read)
			if ablob == None:
				print "Error. Could not read blob measurement info from " + fname
				self.reset()
				return
			else:
				self.blob.append(ablob)

		aworld = FFEA_world_measurement(world_fname, self.num_blobs, num_frames_to_read)
		if aworld == None:
			print "Error. Could not read world measurement info from " + fname
			self.reset()
			return
		else:
			self.world = aworld

		self.num_frames = self.blob[0].num_frames

	def reset(self):

		self.num_blobs = 0
		self.blob = []
		self.world = None

	def plot_energies(self, kT = 1, step_length = 1, num_nodes = None):

		# Get a container and work out the number of columns and rows needed
		plt.figure(1)
		num_columns = 2
		num_rows = self.num_blobs

		# Plot everything
		for i in range(self.num_blobs):

			#
			# Kinetic
			#

			# Measured
			subplt = plt.subplot(num_rows, num_columns, 2 * i)
			subplt.plot(self.blob[i].step * step_length, self.blob[i].ke * 1.0 / kT, label="Trajectory")

			# Theory
			if num_nodes != None and kT != 1:
				subplt.plot(self.blob[i].step * step_length, np.array([3 * num_nodes[i][0] / 2.0 for j in range(self.num_frames)]), "k--", label="Theoretical")

			# Axes
			plt.ylim((0,  max(self.blob[i].ke * 1.0 / kT) * 1.5))

			# Legend
			legend = subplt.legend(loc=1, shadow=True)

			# Titles
			plt.title("Kinetic Energy Trajectory")
			if kT == 1:
				plt.ylabel("Kinetic Energy (J)")
			else:
				plt.ylabel("Kinetic Energy (kT)")

			if step_length == 1:
				plt.xlabel("Step")
			else:
				plt.xlabel("Time (s)")

			#
			# Potential
			#

			# Measured
			subplt = plt.subplot(num_rows, num_columns, 2 * i + 1)
			subplt.plot(self.blob[i].step * step_length, self.blob[i].pe * 1.0 / kT, label="Trajectory")

			# Theory
			if num_nodes != None and kT != 1:
				subplt.plot(self.blob[i].step * step_length, np.array([(3 * num_nodes[i][0] - 6) / 2.0 for j in range(self.num_frames)]), "k--", label="Theoretical")

			# Axes
			plt.ylim((0,  max(self.blob[i].pe * 1.0 / kT) * 1.5))	
			
			# Legend
			legend = subplt.legend(loc=1, shadow=True)

			# Titles
			plt.title("Potential Energy Trajectory")
			if kT == 1:
				plt.ylabel("Potential Energy (J)")
			else:
				plt.ylabel("Potential Energy (kT)")

			if step_length == 1:
				plt.xlabel("Step")
			else:
				plt.xlabel("Time (s)")


			plt.show()

class FFEA_blob_measurement:

	def __init__(self, fname, frames_to_read):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. File " + fname  + " not found."
			return None

		self.num_frames = -1
		just_restarted = 0
		
		print "Reading blob measurement data..."
		for line in fin.readlines():
			if(self.num_frames >= frames_to_read):
				break

			if self.num_frames == -1:
				if line.strip().split("|")[0].strip() != "# step":
					print "Error. Expected '# step | KE | PE | CoM x | CoM y | CoM z | L_x | L_y | L_z | rmsd | vdw_area_%d_surface | vdw_force_%d_surface | vdw_energy_%d_surface' as a header but found " + line
					return None
				else:
					self.num_frames += 1
					continue

			if line.strip() == "#==RESTART==":
				self.num_restarts += 1
				just_restarted = 1
				continue

			if self.num_frames % 100 == 0:
				sys.stdout.write("\tRead " + str(self.num_frames) + " frames")
				if frames_to_read < float("inf"):
					sys.stdout.write(" out of " + str(int(frames_to_read)) + "\n")

			sline = line.split()
			if len(sline) != 13:
				print "Error. Expected 13 columns but found " + str(len(sline))
				self.reset()
				return None

			# Assign stuff
			if just_restarted == 1:
				while(True):
					if self.step[-1] >= int(sline[0]):
						self.num_frames -= 1
						self.ke.pop()
						self.pe.pop()
						self.com.pop()
						self.angmom.pop()
						self.rmsd.pop()
						self.vdw_surfacearea.pop()
						self.vdw_surfaceforce.pop()
						self.vdw_surfaceenergy.pop()

						if self.step[-1] == int(sline[0]):
							self.step.pop()
							break
						else:
							self.step.pop()
					else:
						break

				just_restarted = 0

			self.step.append(int(sline[0]))
			self.ke.append(float(sline[1]))
			self.pe.append(float(sline[2]))
			self.com.append([float(sline[i]) for i in range(3,6)])
			self.angmom.append([float(sline[i]) for i in range(6,9)])
			self.rmsd.append(float(sline[9]))
			self.vdw_surfacearea.append(float(sline[10]))
			self.vdw_surfaceforce.append(float(sline[11]))
			self.vdw_surfaceenergy.append(float(sline[12]))

			self.num_frames += 1

		fin.close()

		print "...done! Read %d frames." % (self.num_frames)

		# Numpy stuff up
		self.step = np.array(self.step)
		self.ke = np.array(self.ke)
		self.pe = np.array(self.pe)
		self.com = np.array(self.com)
		self.angmom = np.array(self.angmom)
		self.rmsd = np.array(self.rmsd)
		self.vdw_surfacearea = np.array(self.vdw_surfacearea)
		self.vdw_surfaceforce = np.array(self.vdw_surfaceforce)
		self.vdw_surfaceenergy = np.array(self.vdw_surfaceenergy)

	def reset(self):

		self.num_frames = 0
		self.num_restarts = 0
		self.step = []
		self.ke = []
		self.pe = []
		self.com = []
		self.angmom = []
		self.rmsd = []
		self.vdw_surfacearea = []
		self.vdw_surfaceforce = []
		self.vdw_surfaceenergy = []

		
	
class FFEA_world_measurement:

	def __init__(self, fname, num_blobs, frames_to_read):

		# Initialise stuff
		self.reset()
		self.num_blobs = num_blobs

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print "Error. File " + fname  + " not found."
			return None

		self.num_frames = -1
		num_columns = self.num_blobs * (self.num_blobs - 1) * 3 / 2 + 1

		print "Reading world measurement data..."
		for line in fin.readlines():

			if(self.num_frames >= frames_to_read):
				break

			if self.num_frames == -1:
			
				# Check correct number of columns
				sline = line.split("|")
				if len(sline) != num_columns:
					print "Error. Expected 3 columns for each possible blob pairing (%d) but got %d columns" % (num_columns, len(sline))
					return None
				else:
					self.num_frames += 1
					continue

			if self.num_frames % 100 == 0:
				sys.stdout.write("\tRead " + str(self.num_frames) + " frames")
				if frames_to_read < float("inf"):
					sys.stdout.write(" out of " + str(int(frames_to_read)) + "\n")

			self.num_frames += 1

		fin.close()

	def reset(self):
		
		self.num_frames = 0
