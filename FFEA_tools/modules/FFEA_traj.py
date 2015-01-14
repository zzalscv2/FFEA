import sys, os
from math import *
from Vectors import *

class FFEA_traj:

	def __init__(self, traj_fname, num_frames_to_read):

		self.num_blobs = 0
		self.blob = []
		self.num_frames = 0
		self.read_traj_from_file(traj_fname, num_frames_to_read)

	def read_traj_from_file(self, traj_fname, num_frames):
		
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
	
		# Get num_nodes for each blob and create traj_blobs
		line = traj.readline().split()
		for i in range(self.num_blobs):
			num_nodes = int(line[4 * i + 3].strip())
			self.blob.append(FFEA_traj_blob(num_nodes))
			print "\tBlob " + str(i) + " has " + str(num_nodes) + " nodes"
	
		# Whitespace
		traj.readline()
		
		# Read trajectory for desired number of frames or all frames
		check = 0
		if num_frames <= 0:
			num_frames = 10000
			
		print("\n\tReading position data")
		for i in range(num_frames):
			
			# Check for asterisk
			if traj.readline().strip() != "*":
				sys.exit("Error. Expected an asterisk to begin frame " + str(i))

			for j in range(self.num_blobs):

				# Check for eof
				line = traj.readline()
				if line == "":
					verify = raw_input("Specified " + str(num_frames) + " to read, but eof reached at " + str(i - 1) + ". Continue (y/n)?:")
					if verify[0] == "y" or verify[0] == "Y":
						print "Continuing"					
						return
					else:
						sys.exit("Error. Trajectory does not contain number of frames specified. Bye!")
		
				# Check for STATIC or not
				line = traj.readline().strip()
				if line == "STATIC":
					continue
			
				# Read frame
				work_frame = FFEA_traj_frame()
				for k in range(self.blob[j].num_nodes):
					line = traj.readline().split()
					work_vector = vector3(float(line[0].strip()), float(line[1].strip()), float(line[2].strip()))
					work_frame.node.append(work_vector)
				
				self.blob[j].frame.append(work_frame)
			
			check += 1
			if num_frames >= 100 and check % (num_frames / 50) == 0:
				print "\t\tRead " + str(check) + " frames"

		print "Done. Read " + str(check) + " frames in total."
		self.num_frames = check
		for i in range(self.num_blobs):
			self.blob[i].num_frames = len(self.blob[i].frame)		

		# Deleting stuff just in case each blob has different number of frames
		print "Deleting extra frames"
		for i in range(self.num_blobs):
			if self.num_frames != self.blob[i].num_frames:
				difference = self.num_frames - self.blob[i].num_frames
				print "\tBlob " + str(i) + " has " + str(difference) + " more frames than others! Deleting..."
				for j in range(difference):
					self.blob[i].frame.pop()

				self.blob[i].num_frames = self.num_frames

		print "Done. All blobs now have " + str(check) + " frames in total."
		return
			
	def calc_centroids(self):

		for ablob in self.blob:
			for aframe in ablob.frame:
				aframe.calc_centroid()
						
	def write_centroids_to_file(self, fname):
		
		fout = open(fname, "w")
		fout.write("FFEA_traj Centroid Trajectories\n\n")
		fout.write("Frame\t")
		for i in range(self.num_blobs):
			fout.write("Blob " + str(i) + " centroid\t")
		
		fout.write("\n")
		for i in range(self.num_frames):
			fout.write(str(i) + "\t")
			for ablob in self.blob:
				fout.write("%5.2e %5.2e %5.2e\t" % (ablob.frame[i].centroid.x, ablob.frame[i].centroid.y, ablob.frame[i].centroid.z))
			fout.write("\n")
		fout.close()

class FFEA_traj_blob:

	def __init__(self, num_nodes):
		
		self.num_nodes = num_nodes
		self.num_frames = 0
		self.frame = []
	
class FFEA_traj_frame:

	def __init__(self):
		
		self.node = []
		self.centroid = vector3(0.0, 0.0, 0.0)

	def calc_centroid(self):
		
		self.centroid.set_pos(0.0, 0.0, 0.0)
		for anode in self.node:
			self.centroid += anode

		self.centroid.scale(1.0/len(self.node))
