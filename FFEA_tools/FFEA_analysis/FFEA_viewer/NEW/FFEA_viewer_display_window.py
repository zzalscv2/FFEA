from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

import threading

#import pyximport; pyximport.install()
from Quaternion import *

import time

import Image

import shutil

import numpy as np

import FFEA_script, FFEA_viewer_blob, FFEA_spring

class FFEA_viewer_display_window:

	def __init__(self, speak_to_control, ffea_fname, energy_thresh=1.0e6):

		self.energy_threshold = energy_thresh
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

		# params object
		self.pms = None

		# camera
		self.orientation = Quaternion()
		self.z = 1

		# mouse
		self.last_x = -1
		self.last_y = -1

		# frames
		self.frame = 0
		self.num_frames = 0

		# list of loaded blobs
		self.blob = []

		# List of global objects
		self.springs = None
		self.animate = False
		self.speed = 1
		self.pause_loading = False
		self.pausing = False

		# Selected stuff
		self.selected_index = 0
		self.selected_blob = 0
		self.selected_conformation = 0

		self.display_flags = {	'show_mesh': 0,
					'show_solid': 1,
					#'show_flat': 0,
					'show_material': 0,
					'show_vdw_only': 0,
					'show_node_numbers': 0,
					'show_pinned_nodes': 1,
					#'hide_frozen': 0,
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


		self.offset_x = 0
		self.offset_y = 0
		self.offset_z = 0

		self.box = False
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

		# Get script as an object
		script = FFEA_script.FFEA_script(self.ffea_fname)
		self.pms = script.params
		
		# Let control window know about num_blobs and num_conformations
		self.speak_to_control.send({'num_blobs': self.pms.num_blobs, 'num_conformations': self.pms.num_conformations})

		# Build box (if possible)
		if(self.pms.es_N_x != 0 and self.pms.es_N_y != 0 and self.pms.es_N_z != 0):
			self.box = True
			self.box_x = (1.0/self.pms.kappa) * self.pms.es_h * self.pms.es_N_x
			self.box_y = (1.0/self.pms.kappa) * self.pms.es_h * self.pms.es_N_y
			self.box_z = (1.0/self.pms.kappa) * self.pms.es_h * self.pms.es_N_z

		# Now, load every blob and conformation into a nice 2D list
		i = 0;
		bi = 0
		ci = 0
		for b in script.blob:
			conf_list = []
			for c in b.conformation:
				conf = FFEA_viewer_blob.FFEA_viewer_blob()
				conf.load(bi, ci, nodefname = c.nodes, surffname = c.surface, topfname = c.topology, matfname = c.material, vdwfname = c.vdw, pinfname = c.pin, stokesfname = c.stokes, bsitesfname = c.bsites, motion_state = c.motion_state)
				conf_list.append(conf)

				# Send some details to the control window
				new_blob_name = ffea_id_string + "#" + str(bi) + ", " + str(ci)
				info_string = "Name:\t" + ffea_id_string + "\nConformation:\t" + str(ci) + "\nNodes:\t" + c.nodes + "\nTopology:\t" + c.topology + "\nSurface:\t" + c.surface + "\nMaterial:\t" + c.material + "\nVdW:\t" + c.vdw + "\nPin:\t" + c.pin + "\nStokes:\t" + c.stokes + "\nBinding Sites:\t" + c.bsites + "\nMotion State:\t" + c.motion_state + "\n"
				add_blob_info = {'name': new_blob_name, 'info': info_string}
				self.speak_to_control.send({'add_blob': add_blob_info})
				i += 1
				ci += 1

			self.blob.append(conf_list)
			ci = 0
			bi += 1
			
		# And now the global things
		self.springs = FFEA_spring.FFEA_springs(script.spring)

		# Put beads here and then write a draw function for them for visualisation :)
		#self.beads = FFEA_beads.FFEA_beads(script.beadfname)

		#
		# Translate and Rescale the system
		#

		# Translate, rotate and rescale the initial node files, to put them on par with the trajectories. Also get the smallest scale in the system
		global_scale = float("inf")
		for b in self.blob:
			bi = self.blob.index(b)
			if script.blob[bi].scale < global_scale:
				global_scale = script.blob[bi].scale

			for c in b:
				c.set_centroid(script.blob[bi].centroid)
				c.apply_rotation(script.blob[bi].rotation)

		for b in self.blob:
			for c in b:
				c.node.pos *= script.blob[bi].scale / global_scale

		# Rescale box
		self.box_x *= global_scale
		self.box_y *= global_scale
		self.box_z *= global_scale


		# Now, move initial structure into simulation into box, if necessary
		if self.pms.move_into_box == 1 and self.pms.calc_vdw == 1:
			world_centroid = np.array([0.0, 0.0, 0.0])
			shift = np.array([0.0, 0.0, 0.0])
			total_num_nodes = 0

			for b in self.blob:

				bcent = b[0].node.calc_centroid()
				world_centroid += bcent * b[0].node.num_nodes
				total_num_nodes += b[0].node.num_nodes

			world_centroid *= 1.0 / total_num_nodes
			shift[0] = self.box_x / 2.0 - world_centroid[0]
			shift[1] = self.box_y / 2.0 - world_centroid[1]
			shift[2] = self.box_z / 2.0 - world_centroid[2]

			for b in self.blob:
				b[0].node.translate(shift)

		# If nodes are STATIC, use this as a frame. If trajectory doesn't exist, use this as a frame
		#if self.pms.trajectory_out_fname == None or not os.path.exists(self.pms.trajectory_out_fname):
		if not self.pms.trajectory_out_fname == None:
			for b in self.blob:
				for c in b:
					if b.index(c) == 0:
						c.set_frame_from_nodes()
					else:
						c.add_empty_frame()

		else:
			for b in self.blob:
				for c in b:
					if b.index(c) == 0 and c.motion_state == "STATIC":
						c.set_frame_from_nodes()

		# Now load trajectory
		#if trajectory_out_fname != None:
		#	self.load_trajectory_thread = threading.Thread(target=self.load_trajectory, args=(trajectory_out_fname,))
		#	self.load_trajectory_thread.start()

		# Hold on calculating dimensions until at least one frame has been calculated from a trajectory, if it exists
		#while(self.num_frames < 1):
		#	if self.pms.trajectory_out_fname == None:
		#		break
		#	else:
		#		pass

		# Reset initial camera (dependent upon structure size)
		dims = self.get_system_dimensions()
		self.dimensions = [dims[i][1] - dims[i][0] for i in range(3)]

		if (self.dimensions[2] > self.dimensions[1]) and (self.dimensions[2] > self.dimensions[0]):
			self.z = 2 * self.dimensions[2]
		elif self.dimensions[0] > self.dimensions[1]:
			self.z = self.dimensions[0] / (2 * np.tan(np.pi / 6.0))
		else:
			self.z = self.dimensions[1] / (2 * np.tan(np.pi / 6.0))


	def load_trajectory(self, trajectory_out_fname,):

		# Firstly, delete the frames loaded from node files originally
		for blob in self.blob:
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
		if self.pms.num_blobs != line:
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
				self.pms.num_conformations[i] = 1

			# Get each num_nodes
			sline = line.split()
			for i in range(self.pms.num_blobs):
				self.blob[i][0].num_nodes = int(sline[4 * i + 3])
		else:
			sline = line.split()[3:]
		
			for i in range(len(sline)):
				if self.pms.num_conformations[i] != int(sline[i]):
					print "Error. 'Number of Conformations' %d specified in trajectory file not consistent with script file." % (i)
					no_errors_loading = False
					traj.close()
					return

			# Get each num_nodes
			for i in range(self.pms.num_blobs):
				sline = traj.readline().split()[2:]
				for j in range(self.pms.num_conformations[i]):
					self.blob[i][j].num_nodes = int(sline[4 * j + 3])
	
		# Final whitespace
		traj.readline()

		# Start Reading Frames
		active_conf = [0 for i in range(self.pms.num_blobs)]
		completed = 0
		trajtype = 1	# New type
		tested = 0
		first_frame = self.pms.num_blobs
		while True:

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

			for i in range(self.pms.num_blobs):

				# Get active blob (may need to change trajectory format in future due to repeated info)
				try:
					filepos = traj.tell()
					active_conf[i] = int(traj.readline().split()[3][0])
					traj.seek(filepos)

				except(IndexError):
					traj.close()
					completed = 1
					break

				# Load actual frame
				j = active_conf[i]

				# if the blob has state STATIC, then there is no node information in the trajectory file (since it is unchanged during the simulation),
				# therefore just read in the word STATIC
				if self.blob[i][j].get_state() == "STATIC":
					traj.readline() # skip "Blob x, Conformation y, step z" line
					traj.readline() # skip "STATIC" line
					continue

				# for DYNAMIC or FROZEN blobs, try to read the node info for this frame
				try:
					self.blob[i][j].load_frame(traj)
					for k in range(self.pms.num_conformations[i]):
						if k != j:
							self.blob[i][k].load_frame(None)

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
							for i in range(self.pms.num_blobs):
								traj.readline()
	
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
		
						for i in range(self.pms.num_blobs):
							traj.readline()

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

	def get_system_dimensions(self):

		dims = [[float("inf"), -1* float("inf")] for i in range(3)]

		for b in self.blob:
			bdims = b[0].get_dimensions()

			for i in range(3):
				if bdims[i][0] < dims[i][0]:
					dims[i][0] = bdims[i][0]
				if bdims[i][1] > dims[i][1]:
					dims[i][1] = bdims[i][1]
		
		return dims

	def death(self):
		glutLeaveMainLoop()

	def change_indices(self, index):

		# Change to control blob index to selected_blob and selected_conformation
		self.selected_index = index
		k = 0
		for i in range(self.pms.num_blobs):
			for j in range(self.pms.num_conformations[i]):
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
		for i in range(self.pms.num_blobs):
			for j in range(self.pms.num_conformations[i]):
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
				self.blob[self.selected_blob][self.selected_conformation].write_vdw(control_stuff['save_vdw'])
			elif "show_blob" in control_stuff.keys():
				index = control_stuff["show_blob"]
				blobi, confi = self.return_changed_indices(index)
				print "Showing blob ", blobi
				for i in range(self.pms.num_conformations[blobi]):
					self.blob[blobi][i].show()

			elif "hide_blob" in control_stuff.keys():
				index = control_stuff["hide_blob"]
				blobi, confi = self.return_changed_indices(index)
				print "Hiding blob ", blobi
				for i in range(self.pms.num_conformations[blobi]):
					self.blob[blobi][i].hide()

			else:	
				if control_stuff['death'] == True:
					self.death()
					return
				self.animate = control_stuff['animate']
				self.speed = control_stuff['speed']
				self.pause_loading = control_stuff['pause_loading']
				self.display_flags = control_stuff['display_flags']

				if self.selected_index != self.display_flags['selected_index']:
					self.change_indices(self.display_flags['selected_index'])
					self.offset_x = 0
					self.offset_y = 0
					self.offset_z = 0

				self.show_box = control_stuff["show_box"]
				if control_stuff['change_frame_to'] != -1:
					self.frame = int(control_stuff['change_frame_to'])
					self.change_indices(self.display_flags['selected_index'])
					if self.blob[self.selected_blob][self.selected_index].hide_blob == False:
						for i in range(self.pms.num_conformations[self.selected_blob]):
							self.blob[self.selected_blob][i].show()

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

		if len(self.blob) > 0:
			for i in range(self.pms.num_conformations[self.selected_blob]):
				cent = self.blob[self.selected_blob][i].traj.frame[self.frame].calc_centroid()
				if cent[0] != None:
					break

			if self.projection == "orthographic":
				self.set_orthographic_projection();

			glEnable(GL_LIGHTING);
			glMatrixMode(GL_MODELVIEW);
			glViewport(0,0,self.width,self.height);
			glLoadIdentity();
			#glTranslated(-self.offset_x, -self.offset_y, -self.offset_z);

			position = [-cent[0] - self.offset_x, -cent[1] - self.offset_y, -cent[2] - self.offset_z - self.z, 1.0];
			glLightfv(GL_LIGHT0, GL_POSITION, position);

			m = self.orientation.construct_matrix();
			m[12] = -self.offset_x
			m[13] = -self.offset_y
			m[14] = -self.z;
			glLoadMatrixd(m);
			glTranslated(-cent[0] - self.offset_x, -cent[1] - self.offset_y, -cent[2] - self.offset_z);

			if self.show_box == 1:
				self.draw_box()

			for i in range(self.pms.num_blobs):

				# First Check if hidden
				hidden = False
				for j in range(self.pms.num_conformations[i]):
					if self.blob[i][j].hide_blob == True:
						hidden = True
						break
				
				# Now draw if not hidden
				if hidden == False:
					for j in range(self.pms.num_conformations[i]):
						self.blob[i][j].draw_frame(self.frame, self.display_flags)

						
		else:
			cent[0] = 0
			cent[1] = 0
			cent[2] = 0

		glutSwapBuffers()		

	def pick(self, mouse_x, mouse_y):
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
		glRenderMode(GL_RENDER)
		if len(self.blob) > 0:
			for i in range(self.pms.num_conformations[self.selected_blob]):
				cent = self.blob[self.selected_blob][i].traj.frame[self.frame].calc_centroid()
				if cent[0] != None:
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
			glTranslated(-cent[0] - self.offset_x, -cent[1] - self.offset_y, -cent[2] - self.offset_z);
			
			self.blob[self.selected_blob][self.selected_conformation].draw_pick_frame(self.frame)

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
		if self.box == False:
			print "No box specified..."
			return

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
					self.blob[self.selected_blob][self.selected_conformation].incr_vdw_face(selected_face)

				if self.display_flags['binding_site_edit_mode'] == 1:
					selected_face = self.pick(x, y)
					self.blob[self.selected_blob][self.selected_conformation].add_face_to_binding_site(selected_face)
		
		if self.mouse_button == 5:
			if self.mouse_state == GLUT_DOWN:
				selected_face = self.pick(x, y)
				self.blob[self.selected_blob][self.selected_conformation].hide_unhide_face(selected_face)

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
		#glOrtho(-self.dimensions[0] * 1000, self.dimensions[0] * 1000, -self.dimensions[1] * 1000, self.dimensions[1] * 1000, -self.dimensions[2] * 1000, self.dimensions[2] * 1000);
		glMatrixMode(GL_MODELVIEW)