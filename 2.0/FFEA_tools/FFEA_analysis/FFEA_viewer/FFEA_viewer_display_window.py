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
		
#		self.set_orthographic_projection()

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
					'selected_blob': 0,
					'show_linear_nodes_only': 0,
					'show_mesh_surf': 0,
					'show_inverted': 0,
					'blob_colour': (1.0, 1.0, 1.0)}

		self.selected_blob = 0

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
				self.active_conformation_index = [0] * self.num_blobs

			elif lvalue == "num_conformations":
				rvalue_split = rvalue[1:-1].split(",")
				for value in rvalue_split:
					self.num_conformations.append(int(value))
					self.total_num_blobs += self.num_conformations[-1]

		self.box_x = (1.0/kappa) * es_h * es_N_x
		self.box_y = (1.0/kappa) * es_h * es_N_y
		self.box_z = (1.0/kappa) * es_h * es_N_z
		if es_N_z != 0:
			scaling_factor = self.z / self.box_z
			scaling_factor_set = 1
		else:
			scaling_factor = 1.0
			scaling_factor_set = 0

		self.box_x *= scaling_factor
		self.box_y *= scaling_factor
		self.box_z *= scaling_factor

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
				blob_mat = []
				blob_stokes = []
				blob_centroid_pos = None
				scale = 1.0

				while True:
					line = ffea_in.readline().strip()[1:-1]
					if line == "/blob":
						break
					elif line == "conformation":
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
							else:
								sys.exit("In " + self.ffea_fname + ", " + rvalue + " is an unexpected rvalue\n")
					
					elif "=" not in line:
						continue
					else:
						lvalue, rvalue = line.split("=")
						lvalue = lvalue.strip()
						rvalue = rvalue.strip()
						if lvalue == "centroid_pos":
							rsplit = rvalue[1:-1].split(",")
							blob_centroid_pos = [float(val) for val in rsplit]
							
						if lvalue == "scale":
							if scaling_factor_set == 0:
								scaling_factor = 1.0 / float(rvalue)
								scaling_factor_set = 1

							scale = scaling_factor
							
						elif lvalue == "initial_conformation":
							self.active_conformation_index[blob_index] = int(rvalue)

				# Load all of the conformations
				for i in range(self.num_conformations[blob_index]):
					print "\nLoading blob " + str(blob_index) + ", conformation " + str(i)
					new_blob = Blob.Blob(energy_thresh=self.energy_threshold)
					new_blob.load(blob_number, blob_index, conformation_index, blob_nodes[i], blob_top[i], blob_surface[i], blob_vdw[i], scale, blob_motion_state[i], blob_pin[i], blob_centroid_pos)

					self.blob_list.append(new_blob)
					new_blob_name = ffea_id_string + "#" + str(blob_index) + ", " + str(conformation_index)
					info_string = "Name:\t" + ffea_id_string + "\nNodes:\t" + blob_nodes[i] + "\nTopology:\t" + str(blob_top[i]) + "\nSurface:\t" + blob_surface[i] + "\nVdW:\t" + str(blob_vdw[i]) + "\npin:\t" + str(blob_pin[i]) + "\nState:\t" + blob_motion_state[i] + "\n"
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
				blob_mat = []
				blob_stokes = []
				blob_centroid_pos = None
				scale = 1.0
			else:
				sys.exit("Expected a blob block in " + self.ffea_fname + " after the <system> tag\n")

		# Load frames from the trajectory file
		if trajectory_out_fname != None:
			# if any of the blobs are STATIC, just load their node positions from the node file
			for blob in self.blob_list:
				if blob.get_state() == "STATIC":
					blob.load_nodes_file_as_frame()

			# Start loading frames for each blob from the trajectory file
			self.load_trajectory_thread = threading.Thread(target=self.load_trajectory, args=(trajectory_out_fname,))
			self.load_trajectory_thread.start()
			#self.load_trajectory(trajectory_out_fname)
		# else just use the nodes files (if traj file not given or found)
		else:
			print "WARNING: Trajectory file is missing. Loading positions from node files."
			for blob in self.blob_list:
				blob.set_scale(1.0)
				blob.load_nodes_file_as_frame()
			

	def load_trajectory(self, trajectory_out_fname):

		print "Reading in trajectory file " + trajectory_out_fname
		traj = open(trajectory_out_fname, "r")

		self.num_frames = 0
		no_errors_loading = True
		#traj.readline()
		#traj.readline()

		# get all initial stuff
		for i in range(3):
			traj.readline()

		# get num_blobs
		if self.num_blobs != int(traj.readline().split()[3].strip()):
			no_errors_loading = False
			traj.close()
			return

		# Get each num_nodes
		sline = traj.readline().split()
		blob_id = 0
		for i in range(self.num_blobs):
			num_nodes = int(sline[4 * i + 3]) 
			for j in range(self.num_conformations[i]):
				blob_id = 0
				for k in range(i):
					blob_id += self.num_conformations[k]
				blob_id += j
				self.blob_list[blob_id].num_nodes = num_nodes

		# Final whitespace
		traj.readline()
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

			for blob in self.blob_list:
				# Only give traj data to active blob
				if blob.conformation_index != self.active_conformation_index[blob.blob_index]:
					continue

				# if the blob has state STATIC, then there is no node information in the trajectory file (since it is unchanged during the simulation),
				# therefore just read in the word STATIC
				if blob.get_state() == "STATIC":
					traj.readline() # skip "Blob x, Conformation y, step z" line
					traj.readline() # skip "STATIC" line
					continue

				# for DYNAMIC or FROZEN blobs, try to read the node info for this frame
				try:
					blob.load_frame(traj)

				except(IndexError):
					no_errors_loading == False
					break			

			if no_errors_loading == True:
				self.num_frames += 1
			else:
				break

		traj.close()

	def death(self):
		glutLeaveMainLoop()

	def update(self, i):
		if self.speak_to_control.poll() == True:
			control_stuff = self.speak_to_control.recv()

			if "save_screenshot" in control_stuff.keys():
				print "Received request to save screenshot to the file", control_stuff['save_screenshot']
				self.TGA_screenshot(control_stuff['save_screenshot'])
			elif "save_vdw" in control_stuff.keys():
				print "Received request to save the VdW profile of blob number", self.selected_blob, " to the file", control_stuff['save_vdw']
				self.blob_list[self.selected_blob].write_vdw(control_stuff['save_vdw'])
			elif "show_blob" in control_stuff.keys():
				index = control_stuff["show_blob"]
				print "Showing blob", index
				self.blob_list[index].show()
			elif "hide_blob" in control_stuff.keys():
				index = control_stuff["hide_blob"]
				print "Hiding blob", index
				self.blob_list[index].hide()
			else:
				if control_stuff['death'] == True:
					self.death()
					return
				self.animate = control_stuff['animate']
				self.speed = control_stuff['speed']
				self.pause_loading = control_stuff['pause_loading']
				self.display_flags = control_stuff['display_flags']
				if self.selected_blob != control_stuff['selected_blob']:
					self.selected_blob = control_stuff['selected_blob']
					self.offset_x = 0
					self.offset_y = 0
					self.offset_z = 0

				self.show_box = control_stuff["show_box"]
				if control_stuff['change_frame_to'] != -1:
					self.frame = int(control_stuff['change_frame_to'])
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
			centroid_x, centroid_y, centroid_z = self.blob_list[self.selected_blob].get_centroid(self.frame)
			
			if self.projection == "orthographic":
				self.set_orthographic_projection();

			glEnable(GL_LIGHTING);
			glMatrixMode(GL_MODELVIEW);
			glViewport(0,0,self.width,self.height);
			glLoadIdentity();
			glTranslated(-self.offset_x, -self.offset_y, -self.offset_z);

			position = [-centroid_x - self.offset_x, -centroid_y - self.offset_y, -centroid_z - self.offset_z - self.z * 10, 1.0];
			glLightfv(GL_LIGHT0, GL_POSITION, position);

			m = self.orientation.construct_matrix();
			m[12] = -self.offset_x
			m[13] = -self.offset_y
			m[14] = -self.z;
			glLoadMatrixd(m);
			glTranslated(-centroid_x, -centroid_y, -centroid_z);

			if self.show_box == 1:
				self.draw_box()

			i = 0
			for blob in self.blob_list:
				blob.draw_frame(self.frame, self.display_flags, i)
				i += 1
		else:
			centroid_x = 0
			centroid_y = 0
			centroid_z = 0

		glutSwapBuffers()		

	def pick(self, mouse_x, mouse_y):
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
		glRenderMode(GL_RENDER)
		if len(self.blob_list) > 0:
			centroid_x, centroid_y, centroid_z = self.blob_list[self.selected_blob].get_centroid(self.frame)
			
			glDisable(GL_LIGHTING);
			glDisable(GL_FOG);
			glShadeModel(GL_FLAT)
			glMatrixMode(GL_MODELVIEW);
			glViewport(0,0,self.width,self.height);
			glLoadIdentity();
			glTranslated(-self.offset_x, -self.offset_y, -self.offset_z);

			m = self.orientation.construct_matrix();
			m[14] = -self.z;
			glLoadMatrixd(m);
			glTranslated(-centroid_x, -centroid_y, -centroid_z);
			
			self.blob_list[self.selected_blob].draw_pick_frame(self.frame)

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
					self.blob_list[self.selected_blob].incr_vdw_face(selected_face)

		if self.mouse_button == 5:
			if self.mouse_state == GLUT_DOWN:
				selected_face = self.pick(x, y)
				self.blob_list[self.selected_blob].hide_unhide_face(selected_face)

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