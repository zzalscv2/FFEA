#
# A module for dealing with FFEA trajectory files. Contains the FFEA_trajectory class, sub-classes and analysis methods
#

import sys, os
from Vectors import *
from math import * 

class FFEA_trajectory:
	
	def __init__(self, traj_fname, num_frames):
		# Splitting filename up for future file searching
		traj_dirname = os.path.dirname(traj_fname)
		trajectory_lines = open(traj_fname, "r")
		trajectory_lines.readline() #FFEA trajectory file
		trajectory_lines.readline()
		trajectory_lines.readline() #Initialisation

		# Getting num_blobs
		self.num_blobs = int(trajectory_lines.readline().split()[3])
		self.num_nodes = 0

		# Getting num_nodes for all blobs
		num_nodes = [0] * self.num_blobs
		sline = trajectory_lines.readline().split()
		for i in range(self.num_blobs):
			num_nodes[i] = int(sline[4 * (i + 1) - 1])
			self.num_nodes += num_nodes[i]

		trajectory_lines.readline()

		# Creating required stuff
		self.blob = [FFEA_traj_blob(num_nodes[i]) for i in range(self.num_blobs)]
		
		# Beginning to read trajectory
		frames_read = 0
		print "Reading trajectory..."
		for i in range(num_frames):
			if self.import_frame(trajectory_lines) == -1:
				break

			frames_read += 1
			if frames_read % floor(num_frames/10.0) == 0:
				print "Read " + str(frames_read) + " frames"

		print "done. Read " + str(frames_read) + " frames from file.\n"
		trajectory_lines.close()

		self.num_frames = frames_read
		for i in range(self.num_blobs):
			if len(self.blob[i].frame) == 0:
				node_fname = raw_input("Blob " + str(i) + " is a STATIC blob. Enter node filename to extract single frame or press return to continue without :")
				if node_fname == "\n":
					print "No frames imported for blob " + str(i) + ".\n"
					continue

				try:
					node_file = open(node_fname, "r")

				except(IOError):
					node_fname = traj_dirname + "/" + node_fname
					print "Node file not found. Trying " + str(node_fname) + ".\n"

					try:
						node_file = open(node_fname, "r")

					except(IOError):
						print "Node file not found. No frames imported for blob " + str(i) + ".\n"
						continue


				self.blob[i].load_single_frame_from_nodes_file(node_file)

	def import_frame(self, trajectory_lines):
		trajectory_lines.readline() #*
		for i in range(self.num_blobs):
			if self.blob[i].import_blob_frame(trajectory_lines) == -1:
				return -1
	
	# Broken!		
	def calc_centroid(self, frame_num):
		print "Calculating total frame centroid..."
		centroid = vector3(0.0, 0.0, 0.0)
		for i in range(self.num_blobs):
			centroid += self.blob[i].frame[frame_num].calc_centroid()
		
		print self.num_nodes
		centroid.scale(1.0/self.num_nodes)
		print "done.\nTotal centroid at frame " + str(frame_num) + " = (" + str(centroid.x) + ", " + str(centroid.y) + ", " + str(centroid.z) + ")\n"
		return centroid
		
class FFEA_traj_blob:
	
	def __init__(self, num_nodes):
		self.num_nodes = num_nodes
		self.num_frames = 0
		self.sub_blob = []
		self.frame = []

	def import_blob_frame(self, trajectory_lines):
		trajectory_lines.readline() # Blob i, step n	
		line = trajectory_lines.readline().strip() #STATIC/DYNAMIC/FROZEN
		if line == "" or line == " " or line == "\n":
			return -1

		if line != "DYNAMIC":
			return 0

		node = []
		trajectory_lines.readline() #num_nodes, delete later
		for i in range(self.num_nodes):
			sline = trajectory_lines.readline().split()
			node.append(vector3(float(sline[0]), float(sline[1]), float(sline[2])))
		self.frame.append(FFEA_traj_blob_frame(node))
		self.num_frames += 1

	def load_single_frame_from_nodes_file(self, node_lines):
		node_lines.readline() # ffea node file
		self.num_nodes = int(node_lines.readline().split()[1])
		num_surface_nodes = int(node_lines.readline().split()[1])
		num_interior_nodes = int(node_lines.readline().split()[1])
		node_lines.readline() # surface nodes:

		node = []
		for i in range(num_surface_nodes):
			sline = node_lines.readline().split()
			node.append(vector3(float(sline[0]), float(sline[1]), float(sline[2])))
	
		node_lines.readline() #interior nodes:
		for i in range(num_interior_nodes):
			sline = node_lines.readline().split()
			node.append(vector3(float(sline[0]), float(sline[1]), float(sline[2])))
		self.frame.append(FFEA_traj_blob_frame(node))

	def scale_frames(self, scale_factor):
		for i in range(len(self.frame)):
			for j in range(len(self.frame[i].node)):
				self.frame[i].node[j].scale(scale_factor)

	# Defining vectors related to blob
	def define_vector_through_nodes(self, frame_num, n1, n2):
		if n1 == n2:
			print "Cannot define vector separating nodes as nodes are equal.\n"
			sys.exit()
		print "Calculating vector separation of nodes " + str(n1) + " and " + str(n2) + " at frame " + str(frame_num) + "..."
		vector_separation = self.frame[frame_num].node[n2] - self.frame[frame_num].node[n1]
		print "done. Vector separation = (" + str(vector_separation.x) + ", " + str(vector_separation.y) + ", " + str(vector_separation.z) + ")\n" 
		return vector_separation
	
	def define_vector_normal_to_plane(self, frame_num, n1, n2, n3):
		print "Calculating unit vector normal to plane defined by nodes " + str(n1) + " and " + str(n2) + " at frame " + str(frame_num) + "..."
		v1 = self.frame[frame_num].node[n2] - self.frame[frame_num].node[n1]
		v2 = self.frame[frame_num].node[n3] - self.frame[frame_num].node[n1]
		normal_vector = vec3_cross_prod(v1, v2)
		if normal_vector.mag() == 0.0:
			print "Plane incorrectly defined. Vector separation of each pair of nodes are not linearly independant. Likely that two input nodes are the same.\n"
			sys.exit()
		normal_vector.normalise()
		print "done. Normal vector = (" + str(normal_vector.x) + ", " + str(normal_vector.y) + ", " + str(normal_vector.z) + ")\n" 
		return normal_vector


	# Sub blob stuff
	def define_sub_blob(self, list_of_node_indices):
		self.sub_blob.append(list_of_node_indices)
		print "Sub blob " + str(len(self.sub_blob)) + " defined. " + str(len(list_of_node_indices)) + " nodes in size\n"

	def define_sub_blob_radially(self, frame_num, origin, min_radius, max_radius):
		print "Calculating Sub blob; list of nodes in the area between " + str(min_radius) + " and " + str(max_radius) + " radially away from origin (" + str(origin.x) + ", " + str(origin.y) + ", " + str(origin.z) + ")..."
		list_of_node_indices = []
		for i in range(self.num_nodes):
			node_radius = (self.frame[frame_num].node[i] - origin).mag()
			if node_radius > min_radius and node_radius < max_radius:
				list_of_node_indices.append(i)

		self.sub_blob.append(list_of_node_indices)
		print "Sub blob " + str(len(self.sub_blob)) + " defined. " + str(len(list_of_node_indices)) + " nodes in size\n"

	def calc_sub_blob_centroid(self, frame_num, sub_blob_num):
		print "Calculating Sub blob " + str(sub_blob_num) + " centroid, from frame " + str(frame_num) + " ..."
		sub_blob_centroid = self.frame[frame_num].calc_sub_blob_centroid(self.sub_blob[sub_blob_num])
		print "done. Sub blob centroid = (" + str(sub_blob_centroid.x) + ", " + str(sub_blob_centroid.y) + ", " + str(sub_blob_centroid.z) + ")\n"
		return sub_blob_centroid

	# Whole blob/sub blob stuff i.e. requiring all frames
	def calc_mean_pos(self):
		list_of_positions = []
		for i in range(0, self.num_nodes):
			list_of_positions.append(vector3(0.0, 0.0, 0.0))
		self.mean_frame = FFEA_traj_blob_frame(list_of_positions)
		check = 0		
		for aframe in self.frame:
			check += 1
			if check % 1 == 0:
				print "Num frames iterated = " + str(check)
 			check2 = 0
			for anode in aframe.node:
				check2 += 1
				if check2 % 1 == 0:
					print "Num nodes iterated = " + str(check2)
				self.mean_frame.node[aframe.node.index(anode)] += vec3_scale(anode, 1.0/self.num_frames)
	
	def calc_centroid_moment_about_x0(self, mom_num):
		moment = 0.0
		check = 0
		for aframe in self.frame:
			acentroid = aframe.calc_centroid()
			if self.frame.index(aframe) == 0:
				x0 = acentroid
			acentroid -= x0
			moment += pow(acentroid.mag(), mom_num)

			check += 1
			if floor(check * 100/self.num_frames) % 10  == 0:
				print str(check * 100/self.num_frames) + "% completed"

		print "Moment " + str(mom_num) + " = " + str(moment/self.num_frames)
		return moment/self.num_frames
	
	def calc_sub_blob_centroid_moment_about_x0(self, mom_num, sub_blob_num):
		moment = 0.0
		check = 0	
		for aframe in self.frame:
			acentroid = aframe.calc_sub_blob_centroid(self.sub_blob[sub_blob_num])
			if self.frame.index(aframe) == 0:
				x0 = acentroid
			acentroid -= x0
			moment += pow(acentroid.mag(), mom_num)

			check += 1
			if floor(check * 100/self.num_frames) % 10  == 0:
				print str(check * 100/self.num_frames) + "% completed"

		print "Moment " + str(mom_num) + " = " + str(moment/self.num_frames)
		return moment/self.num_frames

	def calc_progressive_centroid_moment_about_x0(self, mom_num):
		moment = 0.0
		check = 0
		outfile = open("progressive_mean.out", "w")
		outfile.write("Percent Complete\tMoment " + str(mom_num) + "\n\n")
		for aframe in self.frame:
			acentroid = aframe.calc_centroid()
			if self.frame.index(aframe) == 0:
				x0 = acentroid
			acentroid -= x0
			moment += pow(acentroid.mag(), mom_num)
			check += 1
				
			if floor(check * 100/self.num_frames) % 5 == 0:
				outfile.write(str(check * 100.0/self.num_frames) + "\t\t\t" + str(moment/self.num_frames) + "\n")

			if floor(check * 100/self.num_frames) % 10 == 0:
				print str(check * 100/self.num_frames) + "% completed"

		outfile.close()
		print "Moment " + str(mom_num) + " = " + str(moment/self.num_frames)
		return moment/self.num_frames

class FFEA_traj_blob_frame:
	
	def __init__(self, list_of_positions):
		self.node = list_of_positions

	def calc_centroid(self):
		centroid = vector3(0.0, 0.0, 0.0)
		for anode in self.node:
			centroid += anode
		
		self.centroid = vector3(centroid.x, centroid.y, centroid.z)
		self.centroid.scale(1.0/len(self.node))
		return self.centroid

	def calc_sub_blob_centroid(self, sub_blob_nodes):
		centroid = vector3(0.0, 0.0, 0.0)
		for nodeindex in sub_blob_nodes:
			centroid += self.node[nodeindex]

		centroid.scale(1.0/len(sub_blob_nodes))
		return centroid


# Extra modules
def count_total_frames(traj_fname):
	trajectory_lines = open(traj_fname, "r")
	num_frames = 0

	while True:
		line = trajectory_lines.readline()
		if line == "":
			break

		if line == "*\n":
			num_frames += 1
		
		if num_frames%500 == 0 and num_frames != 0:
			print "Number of Frames read = " + str(num_frames)

	print "Total Number of Frames = " + str(num_frames)
	return num_frames
