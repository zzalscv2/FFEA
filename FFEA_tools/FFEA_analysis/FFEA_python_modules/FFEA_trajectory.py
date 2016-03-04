import numpy as np
from math import isnan
import sys, os
import FFEA_topology
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

class FFEA_trajectory:

	def __init__(self, fname, num_frames_to_read = float("inf"), frame_rate = 1, load_all = 1):

		# Initialise stuff
		self.reset()
		np.seterr(all='raise')
		self.fname = fname
		if self.fname == None:
			return

		self.num_frames_to_read = num_frames_to_read
		self.frame_rate = frame_rate
		
		# Load header information only!
		self.load_header()

		# Default behaviour is to load everything
		self.num_frames_read = 0
		self.num_frames_skipped = 0

		if load_all == 1:
			while self.num_frames < num_frames_to_read:
				print self.num_frames
				if self.load_frame() == 1:
					break

	def load_header(self):

		# Start reading
		try:
			self.traj = open(self.fname, "r")
		
		except(IOError):
			print("Error. File " + self.fname  + " not found.")
			raise IOError

		# Header
		if self.traj.readline().rstrip() != "FFEA_trajectory_file":
			print("Error. Expected to read 'FFEA_trajectory_file'. This may not be an ffea traj file")
			raise IOError

		self.traj.readline()
		if self.traj.readline().rstrip() != "Initialisation:":
			print("Error. Expected to read 'Initialisation:' to begin the initialisation section.")
			raise IOError

		# Blobs and conformations
		try:
			self.num_blobs = int(self.traj.readline().split()[3])

		except(ValueError, IndexError):
			print("Error. Expected:\nNumber of Blobs %d")
			self.reset()
			self.traj.close()
			return -1

		
		sline = self.traj.readline().split()
		if(sline[0].strip() != "Number"):
			
			# Old trajectory type! Initialise all conformation stuff to 0
			self.type = "OLD"
			try:
				self.num_conformations = [1 for i in range(self.num_blobs)]
				self.num_nodes = [[int(sline[4 * i + 3]) for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]
				self.traj.readline()
			except:
				print("Error. Expected:\nBlob 0 Nodes %d Blob 1 Nodes %d ... Blob %d Nodes %d")
				self.reset()
				self.traj.close()
				return
		else:
			try:
				self.num_conformations = [int(sline[3 + i]) for i in range(self.num_blobs)]
				self.num_nodes = [[0 for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

			except(ValueError, IndexError):
				print("Error. Expected:\nNumber of Conformations %d %d ... %d")
				self.reset()
				self.traj.close()
				return
				
			# Nodes
			for i in range(self.num_blobs):
				try:
					sline = self.traj.readline().split()
					for j in range(self.num_conformations[i]):
						self.num_nodes[i][j] = int(sline[5 + 4 * j])
				except(ValueError, IndexError):
					print("Error. Expected:\nBlob %d:	Conformation %d Nodes %d\n.\n.\n.\nBlob %d:	Conformation %d Nodes %d")	
					self.reset()
					self.traj.close()
					return
			self.traj.readline()

		# Initialise the trajectory object
		self.blob = [[FFEA_traj_blob(self.num_nodes[i][j]) for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

		# Read first asterisk
		if self.traj.readline().rstrip() != "*":
			print("Error. Expected to read '*' to begin the trajectory.")
			self.reset()
			self.traj.close()
			return

	def set_header(self, num_blobs, num_conformations, num_nodes):
		self.num_blobs = num_blobs
		self.num_conformations = num_conformations
		self.num_nodes = num_nodes
		self.blob = [[FFEA_traj_blob(self.num_nodes[i][j]) for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

	def set_single_frame(self, node, surf = None):

		for i in range(self.num_blobs):
			self.blob[i][0].num_nodes = node[i].num_nodes
			self.blob[i][0].frame.append(FFEA_traj_blob_frame(node[i].num_nodes))
			self.blob[i][0].frame[-1].pos = node[i].pos

			if surf != None:
				self.blob[i][0].frame[-1].calc_normals(surf[i])

			for j in range(1,self.num_conformations[i]):
				self.blob[i][j].frame.append(None)

	def scale_last_frame(self, scale):
		for b in self.blob:
			for c in b:
				if c.frame[-1] != None:
					c.frame[-1].pos *= scale

	def load_frame(self, surf = None):

		# Set some parameters
		active_conformation = 0

		# And for safety
		start_frame_pos = self.traj.tell()

		# Begin frame
		if(self.num_frames_read + self.num_frames_skipped >= self.num_frames_to_read):
			#print("Trajectory object only allowed to read %d frames. For more, first increase 'num_frames_to_read'." % (self.num_frames_to_read))
			#return 1
			pass

		else:
			# Check if we are to read this frame or skip it
			if (self.num_frames_skipped + self.num_frames_read) % self.frame_rate != 0:
				num_asterisks = 0
				while(num_asterisks != 2):
					line = self.traj.readline()

					# If traj ends, return pointer to beginning of frame
					if line == "" or line == []:
						self.traj.seek(start_frame_pos)
						return 1

					elif("*" in line):
						num_asterisks += 1

				self.num_frames_skipped += 1
				return 0

			# Read the frame otherwise!
			for i in range(self.num_blobs):
						
				# Read data regarding this frame
				try:
					line = self.traj.readline()
					if line == "" or line == []:
						self.traj.seek(start_frame_pos)
						return 1
					
					sline = line.split()
					if int(sline[1].split(",")[0].strip()) != i:
						raise ValueError
						
					# Set active conformation
					if self.type == "NEW":
						active_conformation = int(sline[3].split(",")[0])

				except(ValueError, IndexError):
					print("Error. Expected 'Blob " + str(i) + ", Conformation " + str(active_conformation) + ", step %d', but got:")
					print(line)
					self.reset()
					self.traj.close()
					return

				# Initialise a frame
				frame = FFEA_traj_blob_frame(self.num_nodes[i][active_conformation])

				# See if anything needs to be read?
				if self.traj.readline().strip() == "STATIC":
					self.blob[i][active_conformation].motion_state = "STATIC"
					continue

				else:
					# num_nodes used to be here, maybe
					if self.type == "OLD":
						last_pos = self.traj.tell()

						# If no 'num_nodes', return to last line
						if len(self.traj.readline().split()) != 1:
							self.traj.seek(last_pos)

				for j in range(self.num_nodes[i][active_conformation]):
					try:
						sline = self.traj.readline().split()
						for k in range(3):						
							frame.pos[j][k] = float(sline[k])

					except(ValueError, IndexError):
						print("Error on frame "+str(num_frames_read+num_frames_skipped-1)+" pos "+str(j)+". Expected '%f %f %f' at the very least.")
						print sline
						self.reset()
						self.traj.close()
						return

				# Set the frame in place. Use None if conformation is not active
				for j in range(self.num_conformations[i]):
					if j == active_conformation:
						self.blob[i][active_conformation].frame.append(frame)
					else:
						self.blob[i][j].frame.append(None)

			# Conformation changes
			try:
				line = self.traj.readline().rstrip()
				if line == "" or line == []:
					self.traj.seek(start_frame_pos)
					return 1

				if self.type == "OLD":
					num_frames_read += 1
					return 0

				elif line != "*":
					print("Error. Expected to read '*' to end the frame and begin the conformation changes section.")
					self.reset()
					self.traj.close()
					return 1

				if self.type == "NEW":

					self.traj.readline()	#'Conformation changes:'
					for i in range(self.num_blobs):
						sline = self.traj.readline().split()
						#active_conformation = int(sline[6])
	
					if self.traj.readline().rstrip() != "*":
						print("Error. Expected to read '*' to end the frame and begin the conformation changes section.")
						self.reset()
						self.traj.close()
						return
			except:
				print("Error. Could not read the conformation changes section.")
				self.reset()
				self.traj.close()
				return

			# Average normals at nodes, (if necessary)
			if (surf != None):

				# Check surf objects are compatible
				for i in range(self.num_blobs):
					if len(surf[i]) != self.num_conformations[i]:
						print "Cannot calculate normals, wrong surface arrays provided..."
						return 1
									
				for b in self.blob:
					bi = self.blob.index(b)
					for c in b:
						ci = b.index(c)
						if c.motion_state == "DYNAMIC" and c.frame[-1] != None:
							c.frame[-1].calc_normals(surf[bi][ci])

			self.num_frames_read += 1
			self.num_frames += 1
			if self.num_frames_read % 100 == 0:
				print("\tRead " + str(self.num_frames_read) + " frames")
				if self.num_frames_to_read < float("inf"):
					print(" out of " + str(int(self.num_frames_to_read)))

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

		conf_index = [0 for i in range(self.num_blobs)]

		for i in range(self.num_frames):
			sys.stdout.write("\r\t%d%% of frames written to file" % ((100 * i) / self.num_frames))
			for blob in self.blob:

				blob_index = self.blob.index(blob)
				fout.write("Blob %d, Conformation %d, step %d\n" % (blob_index, conf_index[blob_index], i))
				fout.write(blob[conf_index[blob_index]].motion_state + "\n")
				if blob[conf_index[blob_index]].motion_state == "DYNAMIC":
					for pos in blob[conf_index[blob_index]].frame[i].pos:
						fout.write("%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e\n" % (pos[0], pos[1], pos[2], 0, 0, 0, 0, 0, 0, 0))
		
			# Conformation changes
				
			fout.write("*\nConformation Changes:\n")
			for j in range(self.num_blobs):
				fout.write("Blob %d: Conformation %d -> " % (j,conf_index[j]))
				if(i == self.num_frames - 1):
					fout.write("Conformation %d\n" % (conf_index[j]))
				else:
					for c in self.blob[j]:
						if c.frame[i + 1] != None:
							conf_index[j] = self.blob[j].index(c)
							break

					fout.write("Conformation %d\n" % (conf_index[j]))

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
		self.traj = None
		self.fname = None
		self.num_frames_to_read = 0
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

	def reset():
		self.frame = []
		self.subblob = []
		self.num_subblobs = 0
		
	def get_num_frames(self):
		return len(self.frame)

	def add_empty_frame(self):
		self.num_nodes = 0;
		self.frame.append == None
	
	def add_frame(self, f):
		self.frame.append(f)
		
	def set_frame_from_nodes(self, node):
		self.num_nodes = node.num_nodes
		self.frame.append(FFEA_traj_blob_frame(self.num_nodes))
		self.frame[0].set_from_nodes(node)

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

		centroid = np.array([[0.0,0.0,0.0] for i in range(len(self.frame))])
		i = -1
		for f in self.frame:
			i += 1
			for n in nodes:
				centroid[i] += f.pos[n]
			centroid[i] *= 1.0 / len(nodes)

		return centroid

	def write_frame_as_nodes(self, fname, frame_index, scale):

		fout = open(fname, "w")
		fout.write("ffea node file\nnum_nodes %d\nnum_surface_nodes %d\nnum_interior_nodes %d\nsurface nodes:\n" % (self.num_nodes, self.num_nodes, 0))
		for i in range(self.num_nodes):
			fout.write("%8.6f %8.6f %8.6f\n" % (self.frame[frame_index].pos[i][0] * scale, self.frame[frame_index].pos[i][1] * scale, self.frame[frame_index].pos[i][2] * scale))
		fout.write("interior nodes:\n")
		fout.close()

	def calc_vector_fluctuations(self, angle_type = 0, subindex = [0, 1], graph = 0, fname = None):

		# First, check num_subblobs is ok
		if self.num_subblobs < 2:
			print("Error. Insufficient number of subblobs available.")
			return None
		
		# Get the needed trajectories
		ctraj = [self.get_centroid_trajectory(subblob_index = i) for i in subindex]

		# Get coordinate system
		axis = np.array([[0.0, 0.0, 0.0] for i in range(3)])
		axis[0] = ctraj[1][0] - ctraj[0][0]
		axis[0] *= 1.0 / np.linalg.norm(axis[0])

		axis[1] = np.array([0.0, 1.0, 0.0])

		axis[2] = np.cross(axis[0], axis[1])
		axis[2] *= 1.0 / np.linalg.norm(axis[2])

		axis[1] = np.cross(axis[2], axis[0])
		axis[1] *= 1.0 / np.linalg.norm(axis[1])

		# Get initial projections (r0xy and r0zx, to measure against) (both should be r0 itself)
		r0 = ctraj[1][0] - ctraj[0][0]
		r0 *= 1.0 / np.linalg.norm(r0)
		r0proj = [r0, r0]

		# Now, get angle trajectories
		angtraj = []
		
		for i in range(len(ctraj[1])):

			# Get the vector separation
			r = ctraj[1][i] - ctraj[0][i]

			# Get it's projections (into the two relevent planes xy and zx)
			rproj = [r - np.dot(r, axis[j])*axis[j] for j in [2,1]]
			rprojnorm = [rp * 1.0 / np.linalg.norm(rp) for rp in rproj]

			# Get angle (in whatever units requires) (include negative angles)
			try:
				ang = np.array([np.arccos(np.dot(rprojnorm[j], r0proj[j])) for j in range(2)])

			except(FloatingPointError):
				
				# Check for nans!
				ang = np.array([0.0, 0.0])
				dp = np.dot(rprojnorm[j], r0proj[j])
				for j in range(2):
					if dp > 1.0:
						ang[j] = 0.0
					elif dp < -1.0:
						ang[j] = 180.0 * (-1 ^ np.random.random_integers(0,1))

			if np.dot(np.cross(rprojnorm[0], r0proj[0]), axis[2]) < 0:
				ang[0] *= -1
			if np.dot(np.cross(rprojnorm[1], r0proj[1]), axis[1]) < 0:
				ang[1] *= -1	
			if(angle_type == 1):
				ang *= 180 / np.pi
			
			angtraj.append(ang)

		# Analyse stuff
		angtraj = np.array(angtraj).transpose()
		mean = np.array([np.mean(traj) for traj in angtraj])
		stdev = np.array([np.std(traj) for traj in angtraj])
		err = np.array([std / np.sqrt(len(angtraj[0])) for std in stdev])

		# Plot graphs if required
		if graph == 1 and fname != None:
			
			base, ext = os.path.splitext(fname)
			for i in range(2):

				# Set figure
				plt.figure(i)

				# Build plots (hist and best fit)
				n, bins, patches = plt.hist(angtraj[i], 50, normed = 1, facecolor='green', alpha=0.75, label="Angular Distribution")
				x = np.linspace(mean[i] - 3 * stdev[i], mean[i] + 3 * stdev[i],50)
				y = mlab.normpdf(x, mean[i], stdev[i])
				l, = plt.plot(x, y, 'r--', linewidth=1, label=r"Best Fit Normal Distribution" + "\n" + r"$  \mu = %5.2f \pm %5.2f$" % (mean[i], err[i]) + "\n" + r"$  \sigma = %5.2f$" % (stdev[i]))

				# Prettify the graphs
				if angle_type == 0:
					plt.xlabel("Angle (radians)")
				else:
					plt.xlabel("Angle (degrees)")

				plt.ylabel("Probability")
				handles, labels = plt.gca().get_legend_handles_labels()
				plt.legend(handles[::-1], labels[::-1], loc=1, fontsize=11, fancybox=True, shadow=True)

				if i == 0:
					plt.title("Cytoplasmic Dynein Stalk\nAngular Fluctuations (Hinge)")
					plt.savefig(base + "_hingefluctuations" + ext)
				else:
					plt.title("Cytoplasmic Dynein Stalk\nAngular Fluctuations (Perpendicular)")
					plt.savefig(base + "_perpfluctuations" + ext)

		data = [mean, err, stdev]
		return data

class FFEA_traj_blob_frame:

	def __init__(self, num_nodes):
		self.pos = np.array([[0.0 for i in range(3)] for j in range(num_nodes)])
		self.normal = None

	def calc_normals(self, surf):

		# Normals are averages at each node, rather than at each face

		# Initialise normal array
		self.normal = np.array([[0.0,0.0,0.0] for i in range(surf.num_surface_nodes)])

		# Add normal from each face to nodes on the face
		for f in surf.face:
			a = self.pos[f.n[2]] - self.pos[f.n[1]]
			b = self.pos[f.n[1]] - self.pos[f.n[0]]
			n = np.cross(a,b)
			n *= 1.0 / np.linalg.norm(n)
			for ni in f.n:
				self.normal[ni] += n

		# Now, renormalise to get average normals
		for n in self.normal:
			if n[0] == 0.0:
				continue

			n *= 1.0 / np.linalg.norm(n)

	def set_from_nodes(self, node):
		self.pos = node.pos
		self.normal = node.normal

	def calc_centroid(self):
		return np.mean(self.pos, axis=0)
		
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
	
