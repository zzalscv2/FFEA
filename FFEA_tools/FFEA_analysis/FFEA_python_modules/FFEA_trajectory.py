import numpy as np
import sys
import FFEA_topology

class FFEA_trajectory:

	def __init__(self, fname, num_frames_to_read = float("inf"), frame_rate = 1):
		
		# Initialise stuff
		self.reset()

		# Start reading
		try:
			fin = open(fname, "r")
		
		except(IOError):
			print("Error. File " + fname  + " not found.")
			return

		# Header
		if fin.readline().rstrip() != "FFEA_trajectory_file":
			print("Error. Expected to read 'FFEA_trajectory_file'. This may not be an ffea traj file")
			return

		fin.readline()
		if fin.readline().rstrip() != "Initialisation:":
			print("Error. Expected to read 'Initialisation:' to begin the initialisation section.")
			return

		# Blobs and conformations
		try:
			self.num_blobs = int(fin.readline().split()[3])

		except(ValueError, IndexError):
			print("Error. Expected:\nNumber of Blobs %d")
			self.reset()
			fin.close()
			return

		
		sline = fin.readline().split()
		if(sline[0].strip() != "Number"):
			
			# Old trajectory type! Initialise all conformation stuff to 0
			self.type = "OLD"
			try:
				self.num_conformations = [1 for i in range(self.num_blobs)]
				self.num_nodes = [[int(sline[4 * i + 3]) for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]
				fin.readline()
			except:
				print("Error. Expected:\nBlob 0 Nodes %d Blob 1 Nodes %d ... Blob %d Nodes %d")
				self.reset()
				fin.close()
				return
		else:
			try:
				self.num_conformations = [int(sline[3 + i]) for i in range(self.num_blobs)]
				self.num_nodes = [[0 for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

			except(ValueError, IndexError):
				print("Error. Expected:\nNumber of Conformations %d %d ... %d")
				self.reset()
				fin.close()
				return
				
			# Nodes
			for i in range(self.num_blobs):
				try:
					sline = fin.readline().split()
					for j in range(self.num_conformations[i]):
						self.num_nodes[i][j] = int(sline[5 + 4 * j])
				except(ValueError, IndexError):
					print("Error. Expected:\nBlob %d:	Conformation %d Nodes %d\n.\n.\n.\nBlob %d:	Conformation %d Nodes %d")	
					self.reset()
					fin.close()
					return
			fin.readline()

		# Initialise the trajectory object
		self.blob = [[FFEA_traj_blob(self.num_nodes[i][j]) for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

		# Begin traj
		if fin.readline().rstrip() != "*":
			print("Error. Expected to read '*' to begin the trajectory.")
			self.reset()
			fin.close()
			return

		# Read frames until everything breaks!
		not_completed = True
		num_frames_read = 0
		num_frames_skipped = 0
		active_conformation = [0 for i in range(self.num_blobs)]

		# First one doesn't count
		if(num_frames_to_read < float("inf")):
			num_frames_to_read += 1
		print("Reading FFEA_trajectory...")
		while(num_frames_read + num_frames_skipped < num_frames_to_read):

			# Check if we are to read this frame
			if (num_frames_skipped + num_frames_read) % frame_rate != 0:
				num_asterisks = 0
				while(num_asterisks != 2):
					line = fin.readline()
					if line == "" or line == []:
						not_completed = False
						break
					elif("*" in line):
						num_asterisks += 1

				num_frames_skipped += 1
				continue

			for i in range(self.num_blobs):
				frame = FFEA_traj_blob_frame(self.num_nodes[i][active_conformation[i]])
				try:
					line = fin.readline()
					if line == "" or line == []:
						not_completed = False
						break
					
					sline = line.split()
					if int(sline[1][0]) != i:
						raise ValueError

					if self.type == "NEW" and int(sline[3][0]) != active_conformation[i]:
						raise ValueError
						

				except(ValueError, IndexError):
					print("Error. Expected 'Blob " + str(i) + ", Conformation " + str(active_conformation[i]) + ", step %d', but got:")
					print(line)
					self.reset()
					fin.close()
					return

				if fin.readline().strip() == "STATIC":
					if(num_frames_read + num_frames_skipped == 0):
						self.blob[i][j].motion_state = "STATIC"

					continue

				else:
					# num_nodes used to be here, maybe
					if self.type == "OLD":
						last_pos = fin.tell()

						# If no 'num_nodes', return to last line
						if len(fin.readline().split()) != 1:
							fin.seek(last_pos)
				for j in range(self.num_nodes[i][active_conformation[i]]):
					try:
						sline = fin.readline().split()
						for k in range(3):						
							frame.pos[j][k] = float(sline[k])

					except(ValueError, IndexError):
						print("Error on frame "+str(num_frames_read+num_frames_skipped-1)+" pos "+str(j)+". Expected '%f %f %f' at the very least.")
						print sline
						self.reset()
						fin.close()
						return

				for j in range(self.num_conformations[i]):
					if j == active_conformation[i]:
						self.blob[i][active_conformation[i]].frame.append(frame)
					else:
						self.blob[i][j].frame.append(None)

			# Maybe we've finished
			if(not_completed):

				# Conformation changes
				try:
					line = fin.readline().rstrip()
					if self.type == "OLD" and (line == "" or line == []):
						print line
						num_frames_read += 1
						break

					elif line != "*":
						print("Error. Expected to read '*' to end the frame and begin the conformation changes section.")
						self.reset()
						fin.close()
						return


					if self.type == "NEW":

						fin.readline()	#'Conformation changes:'
						for i in range(self.num_blobs):
							sline = fin.readline().split()
							active_conformation[i] = int(sline[6])
	
						if fin.readline().rstrip() != "*":
							print("Error. Expected to read '*' to end the frame and begin the conformation changes section.")
							self.reset()
							fin.close()
							return
				except:
					print("Error. Could not read the conformation changes section.")
					self.reset()
					fin.close()
					return

				num_frames_read += 1
				if num_frames_read % 100 == 0:
					print("\tRead " + str(num_frames_read) + " frames")
					if num_frames_to_read < float("inf"):
						print(" out of " + str(int(num_frames_to_read)))
			else:
				break

		print("...done!\nRead " + str(num_frames_read) + " frames.\nSkipped " + str(num_frames_skipped) + " frames.\nTotal frames parsed = " + str(num_frames_read + num_frames_skipped))
		self.num_frames = num_frames_read

	def write_linear_to_file(self, fname, top):

		# Firstly, get a list of linear nodes
		linear_nodes = [[[] for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]
		num_linear_nodes = [[0 for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]
		
		i = -1
		for b in top:
			i += 1
			j = -1
			for c in b:
				j += 1
				linear_nodes[i][j] = list(c.get_linear_nodes())
				num_linear_nodes[i][j] = len(linear_nodes[i][j])

		fout = open(fname, "w")

		# Write header info
		fout.write("FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs %d\n" % (self.num_blobs))
		fout.write("Number of Conformations ")
		for i in range(self.num_blobs):
			fout.write("%d " % (self.num_conformations[i]))
		fout.write("\n")
		for i in range(self.num_blobs):
			fout.write("Blob %d:\t" % (i))
			for j in range(self.num_conformations[i]):
				fout.write("Conformation %d Nodes %d " % (j, num_linear_nodes[i][j]))

			fout.write("\n")
		fout.write("\n")

		# Frames
		fout.write("*\n")
		for i in range(self.num_frames):
			for blob in self.blob:
				blob_index = self.blob.index(blob)
				conformation_index = 0
				fout.write("Blob %d, Conformation %d, step 0\n" % (self.blob.index(blob), 0))
				fout.write(blob[0].motion_state + "\n")
				if blob[0].motion_state == "DYNAMIC":
					j = -1
					for pos in blob[0].frame[i].pos:
						j += 1
						# Only write if linear
						if j in linear_nodes[blob_index][conformation_index]:
							fout.write("%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e\n" % (pos[0], pos[1], pos[2], 0, 0, 0, 0, 0, 0, 0))
		
			# Conformation changes		
			fout.write("*\nConformation Changes:\n")
			for j in range(self.num_blobs):
				fout.write("Blob %d: Conformation %d -> Conformation %d\n" % (j, 0, 0))
			fout.write("*\n")

		fout.close()

	def build_from_pdb(self, pdb, scale = 1.0):

		# Reset entire trajectory
		self.reset()
		
		# Set properties
		self.type = "NEW"
		self.num_blobs = pdb.num_blobs
		self.num_conformations = [1 for i in range(self.num_blobs)]
		self.num_nodes = [[b.num_atoms] for b in pdb.blob]
		self.num_frames = pdb.num_frames
		self.blob = [[FFEA_traj_blob(self.num_nodes[i][0])] for i in range(self.num_blobs)]

		# Fill up the blobs
		for i in range(self.num_frames):
			for j in range(self.num_blobs):
				aframe = FFEA_traj_blob_frame(self.blob[j][0].num_nodes)
				aframe.pos = pdb.blob[j].frame[i].pos * scale
				self.blob[j][0].frame.append(aframe)

	def write_to_file(self, fname):

		print "Writing trajectory to " + fname + "..."
		fout = open(fname, "w")

		# Write header info
		fout.write("FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs %d\n" % (self.num_blobs))
		fout.write("Number of Conformations ")
		for i in range(self.num_blobs):
			fout.write("%d " % (self.num_conformations[i]))
		fout.write("\n")
		for i in range(self.num_blobs):
			fout.write("Blob %d:\t" % (i))
			for j in range(self.num_conformations[i]):
				fout.write("Conformation %d Nodes %d " % (j, self.num_nodes[i][j]))

			fout.write("\n")
		fout.write("\n")

		# Frames
		fout.write("*\n")
		for i in range(self.num_frames):
			sys.stdout.write("\r\t%d%% of frames written to file" % ((100 * i) / self.num_frames))
			for blob in self.blob:
				fout.write("Blob %d, Conformation %d, step 0\n" % (self.blob.index(blob), 0))
				fout.write(blob[0].motion_state + "\n")
				if blob[0].motion_state == "DYNAMIC":
					for pos in blob[0].frame[i].pos:
						fout.write("%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e\n" % (pos[0], pos[1], pos[2], 0, 0, 0, 0, 0, 0, 0))
		
			# Conformation changes		
			fout.write("*\nConformation Changes:\n")
			for j in range(self.num_blobs):
				fout.write("Blob %d: Conformation %d -> Conformation %d\n" % (j, 0, 0))
			fout.write("*\n")
		print("\ndone!")
		fout.close()
		
	def write_frame_to_file(self, fname, frame_index):

		fout = open(fname, "w")

		# Write header info
		fout.write("FFEA_trajectory_file\n\nInitialisation:\nNumber of Blobs %d\n" % (self.num_blobs))
		fout.write("Number of Conformations ")
		for i in range(self.num_blobs):
			fout.write("%d " % (self.num_conformations[i]))
		fout.write("\n")
		for i in range(self.num_blobs):
			fout.write("Blob %d:\t" % (i))
			for j in range(self.num_conformations[i]):
				fout.write("Conformation %d Nodes %d " % (j, self.num_nodes[i][j]))

			fout.write("\n")
		fout.write("\n")

		# Frame
		fout.write("*\n")
		for blob in self.blob:
			fout.write("Blob %d, Conformation %d, step 0\n" % (self.blob.index(blob), 0))
			fout.write(blob[0].motion_state + "\n")
			for pos in blob[0].frame[frame_index].pos:
				fout.write("%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e\n" % (pos[0], pos[1], pos[2], 0, 0, 0, 0, 0, 0, 0))
		
		# Conformation changes		
		fout.write("*\nConformation Changes:\n")
		for i in range(self.num_blobs):
			fout.write("Blob %d: Conformation %d -> Conformation %d\n" % (i, 0, 0))
		fout.write("*")
		fout.close()

	def reset(self):
		self.type = "NEW"
		self.num_blobs = 0
		self.num_conformations = []
		self.num_nodes = []
		self.blob = []
		self.num_frames = 0
  
	def calc_distance_between_subblobs(self, blob_index, conformation_index, subblob_1_id, subblob_2_id):
		traj1 = self.blob[blob_index][conformation_index].get_centroid_trajectory(subblob_1_id)
		traj2 = self.blob[blob_index][conformation_index].get_centroid_trajectory(subblob_2_id)
		dists = np.zeros(len(traj1))
		for frame_num in range(len(traj1)):		
			dist = np.sqrt( np.power(traj1[frame_num][0] - traj2[frame_num][0], 2) + np.power(traj1[frame_num][1] - traj2[frame_num][1], 2) + np.power(traj1[frame_num][2] - traj2[frame_num][2], 2) )
			dists[frame_num] = dist
		return np.array(dists)
  
	def get_lever_angle_trajectory(self, node1, node2, node3, blob_index=0, conformation_index=0):
		lever_angle_trajectory=np.zeros([self.num_frames, 3]) # create empty array
		for i in range(self.num_frames): # for each frame
			frame_angles = np.array(self.blob[blob_index][conformation_index].frame[i].get_lever_angle(node1, node2, node3)) # get angular distribution for frame
			lever_angle_trajectory[i,0] = frame_angles[0] # assign values to elements in array accordingly
			lever_angle_trajectory[i,1] = frame_angles[1]
			lever_angle_trajectory[i,2] = frame_angles[2]
		return lever_angle_trajectory   
   
class FFEA_traj_blob:

	def __init__(self, num_nodes):
		self.num_nodes = num_nodes
		self.frame = []
		self.subblob = []
		self.num_subblobs = 0
		self.motion_state = "DYNAMIC"

	def define_subblob(self, indices):
	
		self.subblob.append(indices)
		self.num_subblobs += 1
		#return self.subblob[-1]
		return self.num_subblobs - 1 #returns index of subblob as variable

	def get_centroid_trajectory(self, subblob_index = -1):

		# Total blob
		if subblob_index == -1:
			nodes = range(self.num_nodes)		
		elif self.num_subblobs <= subblob_index:
				print("Error. Blob only contains ", self.num_subblobs, " subblobs.")
				return None
		else:
			nodes = self.subblob[subblob_index]
			print(nodes)

		centroid = np.array([[0.0,0.0,0.0] for i in range(len(self.frame))])
		i = -1
		for f in self.frame:
			i += 1
			for n in nodes:
				centroid[i] += f.pos[n]
			centroid[i] *= 1.0 / self.num_nodes

		return centroid

class FFEA_traj_blob_frame:

	def __init__(self, num_nodes):
		self.pos = np.array([[0.0 for i in range(3)] for j in range(num_nodes)])

	def get_subblob_frame(self, nodes_list):
		subblob_frame = np.zeros([len(nodes_list), 3])
		for i in range(len(nodes_list)):
			subblob_frame[i] = self.pos[i]
		return subblob_frame
	def get_lever_angle(self, point1_no, point2_no, point3_no):
		# get angle between two lines formed by 2d projections of 3 points
 		# with those lines being drawn from points 1 to 2 and points 2 to 3
		def get_angle_between_three_points(first_point, second_point, third_point):
			line_gradient_1 = get_line_gradient(first_point, second_point)
			line_gradient_2 = get_line_gradient(second_point, third_point)
			angle_rads = get_angle_between_two_gradients(line_gradient_1, line_gradient_2)
			return angle_rads
		def get_angle_between_two_gradients(m1, m2): #in radians
			return np.arctan( (m2-m1)/(1 + (m2*m1))  )
		def get_line_gradient(first_point, second_point):
      			# note: this takes 2D points! slice them with slice_axes first!
			line_gradient = (first_point[0] - second_point[0])/(first_point[1] - second_point[1])
			return line_gradient
		def slice_axes(point, axis1, axis2): #turn a 3d point into a 2d one by removing 1 axis (0 for x, 1 for y, 2 for z)
			point_sliced = np.array([point[axis1], point[axis2] ])
			return point_sliced
		def slice_multiple_axes(points_list, axis1, axis2):
			new_points_list = []			
			for point in points_list:
				new_points_list.append(slice_axes(point, axis1, axis2))
			return new_points_list

		points_list = [self.pos[point1_no], self.pos[point2_no], self.pos[point3_no]]
		points_xy = slice_multiple_axes(points_list, 0, 1)
		points_xz = slice_multiple_axes(points_list, 0, 2)
		points_yz = slice_multiple_axes(points_list, 1, 2)
		xy_angle = get_angle_between_three_points(points_xy[0], points_xy[1], points_xy[2])
		xz_angle = get_angle_between_three_points(points_xz[0], points_xz[1], points_xz[2])
		yz_angle = get_angle_between_three_points(points_yz[0], points_yz[1], points_yz[2])
		return xy_angle, xz_angle, yz_angle
		# x= 0, y=1, z=2

# Faster than loading a whole trajectory
def get_num_frames(fname):
	num_asterisks = 0
	with open(fname, "r") as fin:
		for line in fin:
			if line.strip() == "*":
				num_asterisks += 1
	
	return (num_asterisks - 1) / 2
	
