import sys, os
from math import *
import numpy as np

class FFEA_traj:

	def __init__(self, traj_fname, num_frames_to_read, first_frame, last_frame, frame_rate):

		self.num_blobs = 0
		self.total_num_blobs = 0
		self.num_conformations = []
		self.blob = []
		self.num_frames = 0
		self.traj_fname = traj_fname
		self.read_traj_from_file(traj_fname, num_frames_to_read, first_frame, last_frame, frame_rate)

	def read_traj_from_file(self, traj_fname, num_frames, first_frame, last_frame, frame_rate):
		
		print "Reading FFEA trajectory " + traj_fname
		traj = open(traj_fname, "r")
		
		# Get initial crap
		if(traj.readline().strip() != "FFEA_trajectory_file"):
			sys.exit("Error. Expected 'FFEA_trajectory_file' line. This may not be an FFEA trajectory. Bye!")

		for i in range(2):
			traj.readline()

		# Get num_blobs		
		self.num_blobs = int(traj.readline().split()[3].strip())
		print "\tSpecified num_blobs = " + str(self.num_blobs)
	
		# Get num_conformations
		confs = traj.readline().split()[3:]
		for conf in confs:
			self.num_conformations.append(int(conf))
			self.total_num_blobs += self.num_conformations[-1]

		# Get num_nodes for each blob and conformation and create traj_blobs
		self.blob = [[] for i in range(self.num_blobs)]
		for i in range(self.num_blobs):
			sline = traj.readline().split()[2:]
			for j in range(self.num_conformations[i]):
				num_nodes = num_nodes = int(sline[4 * j + 3].strip())
				self.blob[i].append(FFEA_traj_blob(num_nodes))
				print "\tBlob " + str(i) + ", Conformation " + str(j) + " has " + str(num_nodes) + " nodes"

		# Whitespace
		traj.readline()
		
		# Read trajectory for desired number of frames or all frames
		num_frames_parsed = 0
		if num_frames <= 0:
			num_frames = 10000
			
		print("\n\tReading position data")
		completed = 0
		num_frames_read = 0
		current_conformation = [0 for i in range(self.num_blobs)]
		for i in range(num_frames):

			# Check for asterisk
			if traj.readline().strip() != "*":
				sys.exit("Error. Expected an asterisk to begin frame " + str(i))

			for j in range(self.num_blobs):

				# Give inactive conformations a dead frame
				for k in range(self.num_conformations[j]):
					if k != current_conformation[j]:
						self.blob[j][k].frame.append(FFEA_traj_frame(self.blob[j][k].num_nodes, 0))

				# Now read active frame

				# Check for eof
				line = traj.readline()
				if line == "":
					print "Specified " + str(num_frames) + " to read, but eof reached at " + str(i - 1) + "."
					print "Continuing...\n"	
					completed = 1				
					break
		
				# Check for STATIC or not
				line = traj.readline().strip()
				self.blob[j][current_conformation[j]].motion_state = line
				if line == "STATIC":
					continue
			
				# Read frame
				work_frame = FFEA_traj_frame(self.blob[j][current_conformation[j]].num_nodes, 1)
				for k in range(self.blob[j][current_conformation[j]].num_nodes):
					line = traj.readline().split()
					
					for l in range(3):
						work_frame.node_pos[k][l] = float(line[l].strip())
						work_frame.node_vel[k][l] = float(line[l + 3].strip())


				if i < first_frame:
					continue
				elif i > last_frame:
					break

				if i % frame_rate == 0:
					self.blob[j][current_conformation[j]].frame.append(work_frame)

			# If not yet at eof, update progress
			if completed == 0:

				# Get conformation change details
				if traj.readline().strip() != "*":
					sys.exit("Error. Expected '*' to end trajectory data at frame " + str(i))
				if traj.readline().strip() != "Conformation Changes:":
					sys.exit("Error. Expected 'Conformation Changes:' to begin kinetic section at frame " + str(i))

				for j in range(self.num_blobs):
					current_conformation[j] = int(traj.readline().split()[6])

				# Get num_frames read/parsed details
				if i % frame_rate == 0:
					num_frames_read += 1

				num_frames_parsed += 1
				if num_frames >= 100 and num_frames_parsed % (num_frames / 50) == 0:
					print "\t\tParsed " + str(num_frames_parsed) + " frames. Of these, read " + str(num_frames_read) + " frames"
			else:
				break

		print "Done. Parsed " + str(num_frames_parsed) + " frames and stored " + str(num_frames_read) + " frames in total."
		self.num_frames = num_frames_read
		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
				self.blob[i][j].num_frames = len(self.blob[i][j].frame)		

		# Deleting stuff just in case each blob has different number of frames
		print "Deleting extra frames"
		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
				if self.num_frames != self.blob[i][j].num_frames:
					difference = self.num_frames - self.blob[i][j].num_frames
					print "\tBlob " + str(i) + " Conformation " + str(j) + " has " + str(difference) + " more frames than others! Deleting..."
					for k in range(difference):
						self.blob[i][j].frame.pop()

					self.blob[i][j].num_frames = self.num_frames

		print "Done. All blobs now have " + str(num_frames_read) + " frames in total."
		
		# Calculate basic things
		self.calc_centroids()
		return
			
	def calc_centroids(self):

		for ablob in self.blob:
			for conf in ablob:
				for aframe in conf.frame:
					aframe.calc_centroid()


	def write_traj_to_file(self, fname):

		# Open frame_fname
		fout = open(fname, "w")

		# Initial crap
		fout.write("FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs " + str(self.num_blobs) + "\n")
		for i in range(self.num_blobs):
			fout.write("Blob " + str(i) + " Nodes " + str(self.blob[i].num_nodes) + "\t")
		fout.write("\n\n")

		for i in range(self.num_frames):
			fout.write("*\n")
			
			for j in range(self.num_blobs):
				fout.write("Blob " + str(j) + ", Conformation 0, step " + str(i) + "\n")
				fout.write(self.blob[j].motion_state + "\n")
				if self.blob[j].motion_state == "STATIC":
					continue

				for k in range(self.blob[j].num_nodes):
					pos = self.blob[j].frame[i].node_pos[k]
					vel = self.blob[j].frame[i].node_vel[k]
					fout.write("%6.3e %6.3e %6.3e %6.3e %6.3e %6.3e" % (pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]))
			
					for l in range(4):
						fout.write(" %6.3e" % (0.0))
					fout.write("\n")

		fout.write("*\n")
		fout.close()

	def write_frame_to_file(self, frame_num):
		
		# Get frame_fname
		if len(self.traj_fname.split(".")) == 2:
			frame_fname = self.traj_fname.split(".")[0] + "_frame" + str(frame_num) + "." + self.traj_fname.split(".")[1]	
		else:
			frame_fname = self.traj_fname + + "_frame" + str(frame_num) + ".out"

		# Open frame_fname
		fout = open(frame_fname, "w")
		
		# Initial crap
		fout.write("FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs " + str(self.num_blobs) + "\n")
		fout.write("Number of Conformations")
		for i in range(self.num_blobs):
			fout.write(" " + str(self.num_conformations[i]))
		fout.write("\n")
		for i in range(self.num_blobs):
			fout.write("Blob " + str(i) + ": ")
			for j in range(self.num_conformations[i]):
				fout.write("Conformation " + str(j) + " Nodes " + str(self.blob[i][j].num_nodes) + "\t")
			fout.write("\n")

		fout.write("\n")
		
		current_conformation = -1

		# Write single frame
		fout.write("*\n")
		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
				fout.write("Blob " + str(i) + ", Conformation " + str(j) + ", step " + str(frame_num) + "\n")
				fout.write(self.blob[i][j].motion_state + "\n")

				if self.blob[i][j].motion_state == "STATIC":
					continue

				if self.blob[i][j].frame[frame_num].active == 1:
					print "Blob " + str(i) + " was in conformation " + str(j) + " at frame " + str(frame_num) + ".\n"
					current_conformation = j
					for k in range(self.blob[i][j].num_nodes):
						pos = self.blob[i][j].frame[frame_num].node_pos[k]
						vel = self.blob[i][j].frame[frame_num].node_vel[k]
						fout.write("%6.3e %6.3e %6.3e %6.3e %6.3e %6.3e" % (pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]))
			
						for l in range(4):
							fout.write(" %6.3e" % (0.0))
						fout.write("\n")

		# And conformation changes
		fout.write("*\n")
		if current_conformation != -1:
			fout.write("Conformation Changes:\n")
			for i in range(self.num_blobs):
				fout.write("Blob " + str(i) + ": Conformation " + str(current_conformation) + " -> Conformation " + str(current_conformation) + "\n")
			fout.write("*\n")

		fout.close()
		return frame_fname

class FFEA_traj_blob:

	def __init__(self, num_nodes):
		
		self.num_nodes = num_nodes
		self.num_frames = 0
		self.frame = []
		self.motion_state = "DYNAMIC"

class FFEA_traj_frame:

	def __init__(self, num_nodes, active):
		
		self.active = active
		if self.active == 1:
			self.node_pos = np.array([[0.0 for i in range(3)] for j in range(num_nodes)])
			self.node_vel = np.array([[0.0 for i in range(3)] for j in range(num_nodes)])
			self.centroid = np.array([0.0 for i in range(3)])
		else:
			self.node_pos = None
			self.node_vel = None
			self.centroid = None
			
	def calc_centroid(self):
		
		self.centroid = np.mean(self.node_pos, axis=0)

class FFEA_centroid_traj:
	
	def __init__(self, ffea_traj):

		self.num_blobs = ffea_traj.num_blobs
		self.num_frames = ffea_traj.num_frames
		self.pos = np.array([[[0.0 for i in range(self.num_frames)] for j in range(4)] for k in range(self.num_blobs)])
		for i in range(self.num_blobs):
			for j in range(self.num_frames):
				r = 0.0
				for k in range(3):
					self.pos[i][k][j] = ffea_traj.blob[i].frame[j].centroid[k]
					r += pow(self.pos[i][k][j], 2)

				self.pos[i][3][j] = sqrt(r)

	def write_to_file(self, fname):
		print "Printing FFEA centroids to " + fname
		fout = open(fname, "w")
		fout.write("FFEA_traj Centroid Trajectories\n\n")
		fout.write("Frame\t")
		for i in range(self.num_blobs):
			fout.write("Blob " + str(i) + " centroid (x,y,z,r)\t")
		
		fout.write("\n")
		for i in range(self.num_frames):
			fout.write(str(i) + "\t")
			for j in range(self.num_blobs):
				for k in range(4):
					fout.write("%5.2e " % (self.pos[j][k][i]))
				
				fout.write("\t")
			fout.write("\n")
		fout.close()
#
# Non member functions (for dealing with multiple trajectories etc)
#

# Function to calculate number of frames without actually creating a trajectory
def get_num_frames(traj_fname):

	traj = open(traj_fname, "r")
	num_frames = 0
	while(True):
		line = traj.readline()
		if line == "":
			traj.close()
			break

		elif line.strip() == "*":
			num_frames += 1
	
	return num_frames - 1

# Function to return number of blobs without actually creating a trajectory
def get_num_blobs(traj_fname):

	traj = open(traj_fname, "r")
	traj.readline()
	traj.readline()
	traj.readline()
	
	num_blobs = int(traj.readline().split()[3])
	traj.close()

	return num_blobs

def make_single_blob_traj(traj_fname, blob_num, num_frames):

	# Make new_traj_fname
	if len(traj_fname.split(".")) == 2:
		new_traj_fname = traj_fname.split(".")[0] + "_blob" + str(blob_num) + "." + traj_fname.split(".")[1]
	else:
		new_traj_fname = traj_fname + "_blob" + str(blob_num) + ".out"

	# Make trajectory
	traj = FFEA_traj(traj_fname, num_frames)
	num_nodes = traj.blob[blob_num].num_nodes
	num_frames = traj.blob[blob_num].num_frames
	motion_state = traj.blob[blob_num].motion_state

	# Write to new trajectory
	fout = open(new_traj_fname, "w")
	
	# Initial crap
	fout.write("FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs 1\n")
	fout.write("Blob " + str(blob_num) + " Nodes " + str(num_nodes) + "\n\n")

	# Write frames
	for i in range(num_frames):
		fout.write("*\n")
		fout.write("Blob " + str(blob_num) + " Conformation 0, step " + str(i) + "\n")
		fout.write(motion_state + "\n")
		if motion_state == "STATIC":
			continue
		
		for j in range(num_nodes):
			pos = traj.blob[blob_num].frame[i].node_pos[j]
			vel = traj.blob[blob_num].frame[i].node_vel[j]
			fout.write("%6.3e %6.3e %6.3e %6.3e %6.3e %6.3e" % (pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]))
			
		for j in range(4):
			fout.write(" %6.3e" % (0.0))
		fout.write("\n")
		fout.write("*")

	# Close all
	fout.close()
	traj = None

	return new_traj_fname











