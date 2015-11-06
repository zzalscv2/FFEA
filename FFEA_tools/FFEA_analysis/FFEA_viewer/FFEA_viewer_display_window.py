from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

import threading

#import pyximport; pyximport.install()
from Quaternion import *
import Blob

import time

import Image

import shutil

class FFEA_viewer_display_window():

	def __init__(self, speak_to_control, ffea_fname, num_frames_to_read, energy_thresh=1.0e6):

		self.energy_threshold = energy_thresh
		self.num_frames_to_read = num_frames_to_read
		self.speak_to_control = speak_to_control
		self.ffea_fname = ffea_fname
		self.width = 800
		self.height = 800
		self.init_disp()
		self.init_vars()
		self.load_ffea()
		glutMainLoop()

	def init_disp(self):
		glutInit()
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
		glutInitWindowSize(self.width, self.height)
		glutInitWindowPosition(1,30)
		glutCreateWindow("FFEA Viewer - Display: " + self.ffea_fname)
		glClearColor(1.0,1.0,1.0,1.0)

		glShadeModel(GL_SMOOTH)
		glEnable(GL_CULL_FACE)
		glEnable(GL_DEPTH_TEST)

		density = .002
		fog_colour = [.1, .1, .3]
		glEnable(GL_FOG)
		glFogi(GL_FOG_MODE, GL_EXP2)
		glFogfv(GL_FOG_COLOR, fog_colour)
		glFogf(GL_FOG_DENSITY, density)
		glHint(GL_FOG_HINT, GL_NICEST)
		

		global_ambient = [0.4, 0.4, 0.4, 1.0]
		diffuse = [0.3, 0.3, 0.3, 1.0]
		specular = [0.5, 0.5, 0.5, 1.0]
		specular_reflection = [ 1.0, 1.0, 1.0, 1.0 ]
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient)
		glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse)
		glLightfv(GL_LIGHT0, GL_SPECULAR, specular)
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_reflection)
		glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 128)
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
		glEnable(GL_NORMALIZE)
		
		glEnable(GL_LIGHT0);
		
		shininess = 1.0;
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
		glEnable(GL_COLOR_MATERIAL);

		glutDisplayFunc(self.draw_all)
		glutIdleFunc(self.draw_all)
		glutReshapeFunc(self.reshape_handler)
		glutMouseFunc(self.mouse_handler)
		glutMotionFunc(self.mouse_active);
		glutKeyboardFunc(self.keyboard_handler)
		glutTimerFunc(50, self.update, 0)
		glutCloseFunc(self.close_handler)

	def init_vars(self):
		# camera
		self.orientation = Quaternion()
		self.z = 200

		# mouse
		self.last_x = -1
		self.last_y = -1

		# frames
		self.frame = 0
		self.num_frames = 0

		# list of loaded blobs
		self.blob_list = []

		self.animate = False
		self.speed = 1
		self.pause_loading = False
		self.pausing = False

		self.display_flags = {	'show_mesh': 0,
					'show_solid': 1,
					'show_flat': 0,
					'show_vdw_only': 0,
					'show_node_numbers': 0,
					'show_pinned_nodes': 1,
					'hide_frozen': 0,
					'show_shortest_edge': 0,
					'vdw_edit_mode': 0,
					'binding_site_edit_mode': 0,
					'binding_site_list': [],
					'selected_index': 0,
					'selected_blob': 0,
					'selected_conformation':0,
					'show_linear_nodes_only': 0,
					'show_mesh_surf': 0,
					'show_inverted': 0,
					'blob_colour': (1.0, 1.0, 1.0)}

		self.selected_index = 0
		self.selected_blob = 0
		self.selected_conformation = 0

		self.offset_x = 0
		self.offset_y = 0
		self.offset_z = 0

		self.box_x = -1
		self.box_y = -1
		self.box_z = -1

		self.show_box = 0;

		self.modifying_frame = False

		self.recording = 0
		self.movie_dir = "__temp__FFEA_viewer_movie_dir__"

		self.projection = "perspective"

	def load_ffea(self):
		
		print "Loading ffea file: " + self.ffea_fname
		ffea_path, ffea_id_string = os.path.split(self.ffea_fname)
		if ffea_path == "":
			ffea_path = "."

		ffea_in = open(self.ffea_fname, "r")
		
		# Read required stuff from params block
		trajectory_out_fname = None
		kappa = None
		es_N_x = None
		es_N_y = None
		es_N_z = None
		es_h = None
		self.num_blobs = 0
		self.total_num_blobs = 0
		self.num_conformations = []
		while True:

			# Strip tag wrapping 
			line = ffea_in.readline().strip()[1:-1]

			# Check if reached end
			if line == "/param":
				break

			# Shouldn't be the case, but just in case
			if "=" not in line:
				continue

			# Split line
			lvalue, rvalue = line.split("=")
			lvalue = lvalue.strip()
			rvalue = rvalue.strip()

			if lvalue == "trajectory_out_fname":
				trajectory_out_fname = os.path.expanduser(rvalue)
				if os.path.isfile(trajectory_out_fname) == False:
					print "Warning: ", trajectory_out_fname, "doesn't exist."
					trajectory_out_fname = ffea_path + "/" + os.path.split(trajectory_out_fname)[1]
					print "Trying ", trajectory_out_fname, " instead..."
				if os.path.isfile(trajectory_out_fname) == False:
					print "Error: No such file."
					print "Will load nodes file for each Blob instead, if they can actually be found..."
					trajectory_out_fname = None

			elif lvalue == "kappa":
				kappa = float(rvalue)
			elif lvalue == "es_N_x":
				es_N_x = int(rvalue)
			elif lvalue == "es_N_y":
				es_N_y = int(rvalue)
			elif lvalue == "es_N_z":
				es_N_z = int(rvalue)
			elif lvalue == "es_h":
				es_h = int(rvalue)
			elif lvalue == "num_blobs":
				self.num_blobs = int(rvalue)

			elif lvalue == "num_conformations":

				rvalue_split = rvalue[1:-1].split(",")

				# If no traj, only load first conf
				if trajectory_out_fname == None:
					for value in rvalue_split:
						self.num_conformations.append(1)
				else:
					for value in rvalue_split:
						self.num_conformations.append(int(value))

				self.total_num_blobs += self.num_conformations[-1]

				# Error Check
				if len(self.num_conformations) != self.num_blobs:
					sys.exit("Error. Not enough specified 'num_conformations' to create blob array.")

				# Get blob list array
				self.blob_list = [[None for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

		self.box_x = (1.0/kappa) * es_h * es_N_x
		self.box_y = (1.0/kappa) * es_h * es_N_y
		self.box_z = (1.0/kappa) * es_h * es_N_z
	
		# Let control window know about num_blobs and num_conformations
		self.speak_to_control.send({'num_blobs': self.num_blobs, 'num_conformations': self.num_conformations, 'num_frames': self.num_frames, 'current_frame': -1, 'death': False, 'pausing':self.pausing})

		# Now into system block
		while ffea_in.readline().strip() != "<system>":
			continue

		# Must get all stuff for each blob before loading the conformations
		blob_number = 0
		blob_index = 0
		conformation_index = 0

		while True:
			line = ffea_in.readline().strip()[1:-1]
			if line == "":
				continue
			elif line == "/system":
				break
			elif line == "spring":
				# Ignore
				while True:
					line = ffea_in.readline().strip()[1:-1]
					if line == "/spring":
						break
					else:
						continue
			elif line == "blob":
				blob_nodes = []
				blob_top = []
				blob_motion_state = []
				blob_surface = []
				blob_vdw = []
				blob_pin = []
				blob_binding = []
				blob_mat = []
				blob_stokes = []
				blob_centroid_pos = None
				blob_rotation = None
				scale = 1.0

				while True:
					line = ffea_in.readline().strip()[1:-1]
					if line == "/blob":
						break
					elif line == "conformation":
						binding_set = 0
						while True:
							line = ffea_in.readline().strip()[1:-1]
							if line == "/conformation":
								break
							
							lvalue, rvalue = line.split("=")
							lvalue = lvalue.strip()
							rvalue = rvalue.strip()
							if lvalue == "motion_state":
								blob_motion_state.append(rvalue)
								if rvalue == "STATIC":
									blob_top.append("")
									blob_mat.append("")
									blob_stokes.append("")
									blob_pin.append("")
									blob_binding.append("")
							elif lvalue == "nodes":
								blob_nodes.append(rvalue)
							elif lvalue == "topology":
								blob_top.append(rvalue)
							elif lvalue == "surface":
								blob_surface.append(rvalue)
							elif lvalue == "material":
								blob_mat.append(rvalue)
							elif lvalue == "stokes":
								blob_stokes.append(rvalue)
							elif lvalue == "vdw":
								blob_vdw.append(rvalue)
							elif lvalue == "pin":
								blob_pin.append(rvalue)
							elif lvalue == "binding_sites":
								binding_set = 1
								blob_binding.append(rvalue)
                                                        elif lvalue == "beads":
                                                                continue
							else:
								sys.exit("In " + self.ffea_fname + ", " + rvalue + " is an unexpected rvalue\n")
					
						# Binding not vital
						if binding_set == 0:
							blob_binding.append("")

					elif "=" not in line:
						continue
					else:
						lvalue, rvalue = line.split("=")
						lvalue = lvalue.strip()
						rvalue = rvalue.strip()
						print lvalue, rvalue
						if lvalue == "centroid" or lvalue == "centroid_pos":
							rsplit = rvalue[1:-1].split(",")
							blob_centroid_pos = [float(val) for val in rsplit]
							
						if lvalue == "scale":
							scale = float(rvalue)

						if lvalue == "rotation":
							rsplit = rvalue[1:-1].split(",")
							blob_rotation = [float(val) for val in rsplit]

				# Load all of the conformations
				for i in range(self.num_conformations[blob_index]):

					# If no trajectory, only load first ('active') conformations
					if trajectory_out_fname == None and i != 0:
						break

					print "\nLoading blob " + str(blob_index) + ", conformation " + str(i)
					new_blob = Blob.Blob(energy_thresh=self.energy_threshold)
					new_blob.load(blob_number, blob_index, conformation_index, blob_nodes[i], blob_top[i], blob_surface[i], blob_vdw[i], scale, blob_motion_state[i], blob_pin[i], blob_binding[i], blob_centroid_pos, blob_rotation)

					self.blob_list[blob_index][i] = new_blob
					new_blob_name = ffea_id_string + "#" + str(blob_index) + ", " + str(conformation_index)
					info_string = "Name:\t" + ffea_id_string + "\nConformation:\t" + str(i) + "\nNodes:\t" + blob_nodes[i] + "\nTopology:\t" + str(blob_top[i]) + "\nSurface:\t" + blob_surface[i] + "\nVdW:\t" + str(blob_vdw[i]) + "\npin:\t" + str(blob_pin[i]) + "\nMotion State:\t" + blob_motion_state[i] + "\n"
					add_blob_info = {'name': new_blob_name, 'info': info_string}
					self.speak_to_control.send({'add_blob': add_blob_info})
					blob_number += 1
					conformation_index += 1
				
				blob_index += 1
				conformation_index = 0
				blob_nodes = []
				blob_top = []
				blob_motion_state = []
				blob_surface = []
				blob_vdw = []
				blob_pin = []
				blob_binding = []
				blob_mat = []
				blob_stokes = []
				blob_centroid_pos = None
				blob_rotation = None
				scale = 1.0
			elif line == "interactions":
				while ffea_in.readline().strip() != "</interactions>":
					continue
			else:
				sys.exit("Expected a blob block in " + self.ffea_fname + " after the <system> tag\n")

		# Send binding sites to control
		binding_sites = [[0 for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]
		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
				binding_sites[i][j] = self.blob_list[i][j].num_binding_sites
	
		# Get a global scale
		global_scale = float("inf")
		for blob in self.blob_list:
			if blob[0].scale < global_scale:
				global_scale = blob[0].scale
		
		global_scale = 1.0 / global_scale

		# Rescale box
		self.box_x *= global_scale
		self.box_y *= global_scale
		self.box_z *= global_scale
		
		# Load nodes and shift to appropriate positions
		# Firstly, load all nodes from node files to get a global centroid
		world_centroid = [0.0,0.0,0.0]
		total_num_nodes = 0.0
		for blob in self.blob_list:

			# First conformations only
			blob[0].set_scale(global_scale * blob[0].scale)
			blob[0].load_nodes_file_as_frame()
			x,y,z = blob[0].get_centroid(0)
			world_centroid[0] += x * blob[0].num_nodes
			world_centroid[1] += y * blob[0].num_nodes
			world_centroid[2] += z * blob[0].num_nodes
			total_num_nodes += blob[0].num_nodes
		world_centroid[0] *= 1.0 / total_num_nodes
		world_centroid[1] *= 1.0 / total_num_nodes
		world_centroid[2] *= 1.0 / total_num_nodes

		# Translation to box center
		shift = [0.0,0.0,0.0]
		shift[0] = self.box_x / 2.0 - world_centroid[0]
		shift[1] = self.box_y / 2.0 - world_centroid[1]
		shift[2] = self.box_z / 2.0 - world_centroid[2]

		# Shift to the global centroid if static, else load nodes from traj
		if trajectory_out_fname == None:
			for blob in self.blob_list:
				for c in blob:
					if blob.index(c) == 0:
						c.frames[0].translate(shift)
					c.set_scale(global_scale * blob[0].scale)
		
		else:
			for blob in self.blob_list:
				if blob[0].state == "STATIC":
					for c in blob:
						if blob.index(c) == 0:
							c.frames[0].translate(shift)
				else:
					for c in blob:
						c.set_scale(global_scale)

			self.load_trajectory_thread = threading.Thread(target=self.load_trajectory, args=(trajectory_out_fname,))
			self.load_trajectory_thread.start()

		# Load frames from the trajectory file
		#if trajectory_out_fname != None:
		#	# if any of the blobs are STATIC, just load their node positions from the node file
		#	for blob in self.blob_list:
		#		for conf in blob:
		#			if conf.get_state() == "STATIC":
		#				conf.set_scale(global_scale * conf.scale)
		#				conf.load_nodes_file_as_frame()
		#			else:
		#				conf.set_scale(global_scale)
			# Start loading frames for each blob from the trajectory file
#
#			self.load_trajectory_thread = threading.Thread(target=self.load_trajectory, args=(trajectory_out_fname,))
#			self.load_trajectory_thread.start()
			#self.load_trajectory(trajectory_out_fname)

		# else just use the nodes files (if traj file not given or found)
#		else:
#			print "WARNING: Trajectory file is missing. Loading positions from node files."
#			for blob in self.blob_list:
#				#for conf in blob:
#				blob[0].set_scale(global_scale * blob[0].scale)
#				blob[0].load_nodes_file_as_frame()

	def load_trajectory(self, trajectory_out_fname,):

		# Firstly, delete the frames loaded from node files originally
		for blob in self.blob_list:
			if blob[0].get_state() == "DYNAMIC":
				blob[0].delete_all_frames()

		print "Reading in trajectory file " + trajectory_out_fname
		traj = open(trajectory_out_fname, "r")

		trajtype = 1	# New type
		tested = 0

		self.num_frames = -1
		no_errors_loading = True

		# get all initial stuff
		line = traj.readline().strip()
		if(line != "FFEA_trajectory_file"):
			print "Error. Expected 'FFEA_trajectory_file' , but got " + line
			no_errors_loading = False
			traj.close()
			return

		for i in range(2):
			traj.readline()
			
		# Get num_blobs
		line = int(traj.readline().split()[3].strip())
		if self.num_blobs != line:
			print "Error. 'Number of Blobs' specified in trajectory file not consistent with script file."
			no_errors_loading = False
			traj.close()
			return

		# Get num_conformations
		line = traj.readline()
		if line.split()[0].strip() == "Blob":

			# Old type
			trajtype = 0
			for i in range(len(self.num_conformations)):
				self.num_conformations[i] = 1

			# Get each num_nodes
			sline = line.split()
			for i in range(self.num_blobs):
				self.blob_list[i][0].num_nodes = int(sline[4 * i + 3])
		else:
			sline = line.split()[3:]
		
			for i in range(len(sline)):
				if self.num_conformations[i] != int(sline[i]):
					print "Error. 'Number of Conformations' %d specified in trajectory file not consistent with script file." % (i)
					no_errors_loading = False
					traj.close()
					return

			# Get each num_nodes
			for i in range(self.num_blobs):
				sline = traj.readline().split()[2:]
				for j in range(self.num_conformations[i]):
					self.blob_list[i][j].num_nodes = int(sline[4 * j + 3])
	
		# Final whitespace
		traj.readline()

		# Start Reading Frames
		active_conf = [0 for i in range(self.num_blobs)]
		completed = 0
		trajtype = 1	# New type
		tested = 0
		first_frame = self.num_blobs
		while True:
			if self.num_frames >= self.num_frames_to_read:
				break

			while self.pause_loading == True:
				self.pausing = True
				time.sleep(2)
			self.pausing = False

			# skip asterisk
			line = traj.readline().rstrip()
			if line != "*":
				print "Missing '*' at start of frame", self.num_frames
				print "Instead found '" + line + "'"
				break

			for i in range(self.num_blobs):

				# Load actual frame
				j = active_conf[i]
				
				# if the blob has state STATIC, then there is no node information in the trajectory file (since it is unchanged during the simulation),
				# therefore just read in the word STATIC
				if self.blob_list[i][j].get_state() == "STATIC":
					traj.readline() # skip "Blob x, Conformation y, step z" line
					traj.readline() # skip "STATIC" line
					continue

				# for DYNAMIC or FROZEN blobs, try to read the node info for this frame
				try:
					self.blob_list[i][j].load_frame(traj)
					for k in range(self.num_conformations[i]):
						if k != j:
							self.blob_list[i][k].load_frame(None)

				except(IndexError):
					traj.close()
					completed = 1
					break
			
			if completed == 0:

				# Get conformation data!
				
				if trajtype == 1:
					if tested == 0:
			
						# Testing traj version
						tested = 1
						filepos = traj.tell()
						if traj.readline().strip() != "*":
							print "Error. Expected '*' to end trajectory data."
							no_errors_loading = False
							traj.close()
							return

						if traj.readline().strip() != "Conformation Changes:":
							
							# Old type
							trajtype = 0

							# Set back a line, update and finish
							traj.seek(filepos)
							self.num_frames += 1
						else:
							for i in range(self.num_blobs):
								active_conf[i] = int(traj.readline().split()[6])
	
							if no_errors_loading == True:
								self.num_frames += 1
							else:
								break
					else:
						if traj.readline().strip() != "*":
							print "Error. Expected '*' to end trajectory data."
							no_errors_loading = False
							traj.close()
							return
	
						if traj.readline().strip() != "Conformation Changes:":
							print "Error. Expected 'Conformation Changes:' to begin conformation switching data."
							no_errors_loading = False
							traj.close()
							return
		
						for i in range(self.num_blobs):
							active_conf[i] = int(traj.readline().split()[6])
	
						if no_errors_loading == True:
							self.num_frames += 1
						else:
							break
				else:
					if no_errors_loading == True:
						self.num_frames += 1
					else:
						break
			else:
				break

		traj.close()

	def death(self):
		glutLeaveMainLoop()

	def change_indices(self, index):

		# Change to control blob index to selected_blob and selected_conformation
		self.selected_index = index
		k = 0
		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
				if self.selected_index == k:
					self.selected_blob = i
					self.display_flags['selected_blob'] = i
					self.selected_conformation = j
					self.display_flags['selected_conformation'] = j
					return
				else:
					k += 1
 
	def return_changed_indices(self, index):

		# Change to control blob index to selected_blob and selected_conformation
		self.selected_index = index
		k = 0
		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
				if self.selected_index == k:
					selected_blob = i
					selected_conformation = j
					return selected_blob, selected_conformation
				else:
					k += 1
	def update(self, i):
		self.change_indices(self.selected_index)
		if self.speak_to_control.poll() == True:
			control_stuff = self.speak_to_control.recv()

			if "save_screenshot" in control_stuff.keys():
				print "Received request to save screenshot to the file", control_stuff['save_screenshot']
				self.TGA_screenshot(control_stuff['save_screenshot'])
			elif "save_vdw" in control_stuff.keys():
				print "Received request to save the VdW profile of blob number", self.selected_index, " to the file", control_stuff['save_vdw']
				self.blob_list[self.selected_blob][self.selected_conformation].write_vdw(control_stuff['save_vdw'])
			elif "show_blob" in control_stuff.keys():
				index = control_stuff["show_blob"]
				blobi, confi = self.return_changed_indices(index)
				print "Showing blob ", blobi
				for i in range(self.num_conformations[blobi]):
					self.blob_list[blobi][i].show()

			elif "hide_blob" in control_stuff.keys():
				index = control_stuff["hide_blob"]
				blobi, confi = self.return_changed_indices(index)
				print "Hiding blob ", blobi
				for i in range(self.num_conformations[blobi]):
					self.blob_list[blobi][i].hide()

			else:	
				if control_stuff['death'] == True:
					self.death()
					return
				self.animate = control_stuff['animate']
				self.speed = control_stuff['speed']
				self.pause_loading = control_stuff['pause_loading']
				self.display_flags = control_stuff['display_flags']

				if self.selected_index != control_stuff['selected_index']:
					self.change_indices(control_stuff['selected_index'])
					self.offset_x = 0
					self.offset_y = 0
					self.offset_z = 0

				self.show_box = control_stuff["show_box"]
				if control_stuff['change_frame_to'] != -1:
					self.frame = int(control_stuff['change_frame_to'])
					self.change_indices(control_stuff['selected_index'])
					if self.blob_list[self.selected_blob][self.selected_index].hide_blob == False:
						for i in range(self.num_conformations[self.selected_blob]):
							self.blob_list[self.selected_blob][i].show()

					self.modifying_frame = True

				if self.recording != control_stuff["recording"]:
					self.recording = control_stuff["recording"]
					if self.recording == 1:
						print "Starting recording..."
						self.rec_frame = 0
						if os.path.exists(self.movie_dir):
							shutil.rmtree(self.movie_dir)
						os.makedirs(self.movie_dir)
					else:
						print "Ending recording..."
				if self.projection != control_stuff["projection"]:
					self.projection = control_stuff["projection"]
					print "Changing projection matrix to: " + self.projection
					if self.projection == "orthographic":
						self.set_orthographic_projection()
					elif self.projection == "perspective":
						self.set_perspective_projection()
					else:
						print "Unrecognised projection type: " + self.projection

		if self.animate == True:
			self.change_indices(self.selected_index)
			self.frame += self.speed
			if self.frame > self.num_frames:
				self.frame = 0
			if self.recording == 1:
				frame_fname = self.movie_dir + "/%09d.png" % self.rec_frame
				print frame_fname
				self.TGA_screenshot(frame_fname)
				self.rec_frame += 1

		if self.modifying_frame == False:
			self.speak_to_control.send({'num_frames': self.num_frames, 'current_frame': self.frame, 'death': False, 'pausing': self.pausing})
		else:
			self.speak_to_control.send({'num_frames': self.num_frames, 'current_frame': -1, 'death': False, 'pausing':self.pausing})
			self.modifying_frame = False

		glutTimerFunc(100, self.update, 0)
 
	def draw_all(self):
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)

		self.draw_axes()

		if len(self.blob_list) > 0:
			for i in range(self.num_conformations[self.selected_blob]):
				centroid_x, centroid_y, centroid_z = self.blob_list[self.selected_blob][i].get_centroid(self.frame)

				if centroid_x != None:
					break

			if self.projection == "orthographic":
				self.set_orthographic_projection();

			glEnable(GL_LIGHTING);
			glMatrixMode(GL_MODELVIEW);
			glViewport(0,0,self.width,self.height);
			glLoadIdentity();
			#glTranslated(-self.offset_x, -self.offset_y, -self.offset_z);

			position = [-centroid_x - self.offset_x, -centroid_y - self.offset_y, -centroid_z - self.offset_z - self.z * 10, 1.0];
			#position = [centroid_x, centroid_y, centroid_z - self.z * 10, 1.0];
			glLightfv(GL_LIGHT0, GL_POSITION, position);

			m = self.orientation.construct_matrix();
			m[12] = -self.offset_x
			m[13] = -self.offset_y
			m[14] = -self.z;
			glLoadMatrixd(m);
			glTranslated(-centroid_x, -centroid_y, -centroid_z);
			if self.show_box == 1:
				self.draw_box()

			for i in range(self.num_blobs):

				# First Check if hidden
				hidden = False
				for j in range(self.num_conformations[i]):
					if self.blob_list[i][j].hide_blob == True:
						hidden = True
						break
				
				# Now draw if not hidden
				if hidden == False:
					for j in range(self.num_conformations[i]):
						print i, j
						self.blob_list[i][j].draw_frame(self.frame, self.display_flags)

						
		else:
			centroid_x = 0
			centroid_y = 0
			centroid_z = 0

		glutSwapBuffers()		

	def pick(self, mouse_x, mouse_y):
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
		glRenderMode(GL_RENDER)
		if len(self.blob_list) > 0:
			for i in range(self.num_conformations[self.selected_blob]):
				centroid_x, centroid_y, centroid_z = self.blob_list[self.selected_blob][i].get_centroid(self.frame)
				if centroid_x != None:
					break
			
			glDisable(GL_LIGHTING);
			glDisable(GL_FOG);
			glShadeModel(GL_FLAT)
			glMatrixMode(GL_MODELVIEW);
			glViewport(0,0,self.width,self.height);
			glLoadIdentity();
			glTranslated(-self.offset_x, -self.offset_y, -self.offset_z);

			m = self.orientation.construct_matrix();
			m[12] = -self.offset_x
                        m[13] = -self.offset_y
			m[14] = -self.z;
			glLoadMatrixd(m);
			glTranslated(-centroid_x, -centroid_y, -centroid_z);
			
			self.blob_list[self.selected_blob][self.selected_conformation].draw_pick_frame(self.frame)

		glReadBuffer(GL_BACK)
		face_colour = glReadPixels(mouse_x, self.height - mouse_y, 1, 1, GL_RGB, GL_FLOAT)
		pick_r = int(face_colour[0][0][0] * 255.0)
		pick_g = int(face_colour[0][0][1] * 255.0)
		pick_b = int(face_colour[0][0][2] * 255.0)
		selected_face = (pick_r - 1) + pick_g * 255 + pick_b * 255 * 255

		glEnable(GL_LIGHTING)
		glEnable(GL_FOG)
		glShadeModel(GL_SMOOTH)

		return selected_face

	def draw_axes(self):
		self.set_perspective_projection()
		glDisable(GL_LIGHTING);
		glMatrixMode(GL_MODELVIEW);
		glViewport(0,0,100,100);
		glLoadIdentity();
		m = self.orientation.construct_matrix();
		m[14] = -2;
		glLoadMatrixd(m);
		
		glLineWidth(5.0);
		glBegin(GL_LINES);
		glColor3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glColor3f(1.0, 0.0, 0.0);
		glVertex3f(1.0, 0.0, 0.0);
		
		glColor3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glColor3f(0.0, 1.0, 0.0);
		glVertex3f(0.0, 1.0, 0.0);
		
		glColor3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, 0.0);
		glColor3f(0.0, 0.0, 1.0);
		glVertex3f(0.0, 0.0, 1.0);
		glEnd();
		glLineWidth(1.0);
	
		glColor3f(1.0, .7, .7);
		glRasterPos3f(1.0, 0.0, 0.0);
		glutBitmapString(GLUT_BITMAP_HELVETICA_18, "X");
		
		glColor3f(.7, 1.0, .7);
		glRasterPos3f(0.0, 1.0, 0.0);
		glutBitmapString(GLUT_BITMAP_HELVETICA_18, "Y");
		
		glColor3f(.7, .7, 1.0);
		glRasterPos3f(0.0, 0.0, 1.0);
		glutBitmapString(GLUT_BITMAP_HELVETICA_18, "Z");

	def draw_box(self):
		glDisable(GL_LIGHTING)
		glDisable(GL_FOG)
		glBegin(GL_LINE_STRIP)
		glColor3f(0.3, 1.0, 0.3)
		glVertex3f(0.0, 0.0, 0.0)
		glVertex3f(self.box_x, 0.0, 0.0)
		glVertex3f(self.box_x, self.box_y, 0.0)
		glVertex3f(0.0, self.box_y, 0.0)
		glVertex3f(0.0, 0.0, 0.0)
		glVertex3f(0.0, 0.0, self.box_z)
		glVertex3f(self.box_x, 0.0, self.box_z)
		glVertex3f(self.box_x, self.box_y, self.box_z)
		glVertex3f(0.0, self.box_y, self.box_z)
		glVertex3f(0.0, 0.0, self.box_z)
		glEnd()

		glBegin(GL_LINES)
		glVertex3f(self.box_x, 0.0, 0.0)
		glVertex3f(self.box_x, 0.0, self.box_z)
		glVertex3f(self.box_x, self.box_y, 0.0)
		glVertex3f(self.box_x, self.box_y, self.box_z)
		glVertex3f(0.0, self.box_y, 0.0)
		glVertex3f(0.0, self.box_y, self.box_z)
		glEnd()
		glEnable(GL_LIGHTING)
		glEnable(GL_FOG)

	def mouse_handler(self, button, state, x, y):
		self.mouse_state = state;
		self.mouse_button = button;

		if self.mouse_button == 2:
			if self.mouse_state == GLUT_DOWN:
				if self.display_flags['vdw_edit_mode'] == 1:
					selected_face = self.pick(x, y)
					self.blob_list[self.selected_blob][self.selected_conformation].incr_vdw_face(selected_face)

				if self.display_flags['binding_site_edit_mode'] == 1:
					selected_face = self.pick(x, y)
					self.blob_list[self.selected_blob][self.selected_conformation].add_face_to_binding_site(selected_face)
		
		if self.mouse_button == 5:
			if self.mouse_state == GLUT_DOWN:
				selected_face = self.pick(x, y)
				self.blob_list[self.selected_blob][self.selected_conformation].hide_unhide_face(selected_face)

		if self.mouse_state != GLUT_DOWN:
			self.last_x = -1;
			self.last_y = -1;
		
		if self.mouse_button == 3:
			self.z -= self.z*.05;
		elif self.mouse_button == 4:
			self.z += self.z*.05;

	def mouse_active(self, x, y):
		if self.mouse_button == 0:
			if self.last_x == -1:
				self.last_x = x;
				self.last_y = y;
			else:
				self.orientation.rotate(2 * float(x - self.last_x)/self.width * 3.0, 0, 1, 0);
				self.orientation.rotate(2 * float(y - self.last_y)/self.height * 3.0, 1, 0, 0);
				self.last_x = x;
				self.last_y = y;
		elif self.mouse_button == 2:
			if self.last_x == -1:
				self.last_x = x;
				self.last_y = y;
			else:
				self.offset_x -= (x - self.last_x)
				self.offset_y += (y - self.last_y)
				self.last_x = x;
				self.last_y = y;


	def reshape_handler(self, new_width, new_height):
		self.width = new_width
		self.height = new_height
		
	def keyboard_handler(self, key, x, y):
		if key == '\x1b':
			print "Display received ESC key event."
			self.close_handler()		

	def close_handler(self):
		self.speak_to_control.send({'num_frames': self.num_frames, 'current_frame': self.frame, 'death': True})
		self.death()

	def TGA_screenshot(self, fname):
		screenshot = glReadPixels( 0,0, self.width, self.height, GL_RGBA, GL_UNSIGNED_BYTE)
		im = Image.frombuffer("RGBA", (self.width,self.height), screenshot, "raw", "RGBA", 0, 0)
		im.save(fname)

	def set_perspective_projection(self):
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-1, 1, -1, 1, 1, 10000);
		glMatrixMode(GL_MODELVIEW)

	def set_orthographic_projection(self):
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-self.z, self.z, -self.z, self.z, -100, 10000);
		glMatrixMode(GL_MODELVIEW)
