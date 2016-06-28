import sys, os, time
import numpy as np

from pymol import cmd
from pymol.callback import Callback

from Tkinter import *
import tkFileDialog
import tkMessageBox
import tkColorChooser

import Blob
import threading

# PyMOL stuff:
import subprocess, traceback, Pmw

# Temporary solution to take comments out:
import StringIO

# from multiprocessing import Process, Pipe

# do Ben's springs:
import FFEA_springs

# PyMOL stuff:
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain

# FFEA stuff
import FFEA_script
import FFEA_trajectory

def __init__(self):
  """ 
  Init PyMOL, by adding FFEA stuff to the GUI under Plugins
  """ 
  
  self.menuBar.addmenuitem('Plugin', 'command', 
                           'FFEA Loader', label = 'FFEA Loader...', 
                           command = lambda s=self: FFEA_viewer_control_window(s))


class FFEA_viewer_control_window:
  # # # # # # # # # # # # # # # # # # # # # #
  # # Main control window description # # # #
  # # # # # # # # # # # # # # # # # # # # # # 
  def __init__(self, app):
     self.parent = app.root

     self.root = Tk()

     self.root.geometry("700x200")

     self.root.title("FFEA")

     top_frame = Frame(self.root)
     top_frame.pack()

     menubar = Menu(top_frame)

     filemenu = Menu(menubar, tearoff=0)
     filemenu.add_command(label="Load 'ffea' file", command=self.choose_ffea_file_to_load)
     menubar.add_cascade(label="File", menu=filemenu)
     self.root.config(menu=menubar)

     # PLUGIN (separated into mutually exclusive sets. Devs take note!)
     self.show_solid = IntVar()

     self.show_mesh = IntVar()

     self.show_numbers = IntVar()

     self.show_pinned = IntVar()

     self.show_shortest_edge = IntVar()

     self.show_box = IntVar()

     self.show_springs = IntVar()

     self.do_load_trajectory = IntVar()


     self.init_vars()
 
     # # Display flags frame
     display_flags_frame = Frame(self.root, relief=SUNKEN, bd=1)
     display_flags_frame.pack(anchor=CENTER, expand=True)

     #
     # set val for radiobuttons
     #


     # show solid:
     check_button_show_solid = Radiobutton(display_flags_frame, text="Plain Solid", variable=self.show_solid, value=1, command=lambda:self.update_display_flags("show_solid", val=1))
     check_button_show_solid.grid(row=0, column=0)
     check_button_show_solid.select() # that has to match with the default value 1! 

     # show material: 
     check_button_show_material = Radiobutton(display_flags_frame, text="Material", variable=self.show_solid, value=2, command=lambda:self.update_display_flags("show_solid", val=2))
     check_button_show_material.grid(row=0, column=1)

     # show no solid:
     check_button_show_no_solid = Radiobutton(display_flags_frame, text="No Solid", variable=self.show_solid, value=0, command=lambda:self.update_display_flags("show_solid", val=0))
     check_button_show_no_solid.grid(row=0, column=2)

     # show surface mesh:
     check_button_show_mesh_surf = Radiobutton(display_flags_frame, text="Surface Mesh", variable=self.show_mesh, value=2, command=lambda:self.update_display_flags("show_mesh", val=2))
     check_button_show_mesh_surf.grid(row=1, column=0)

     # show whole mesh:
     check_button_show_mesh = Radiobutton(display_flags_frame, text="Whole Mesh", variable=self.show_mesh, value=1, command=lambda:self.update_display_flags("show_mesh", val=1))
     check_button_show_mesh.grid(row=1, column=1)

     # show no mesh:
     check_button_show_no_mesh = Radiobutton(display_flags_frame, text="No Mesh", variable=self.show_mesh, value=0, command=lambda:self.update_display_flags("show_mesh", val=0))
     check_button_show_no_mesh.grid(row=1, column=2)
     check_button_show_no_mesh.select()

     # show node numbers: 
     check_button_show_node_numbers = Radiobutton(display_flags_frame, text="Node Indices", variable=self.show_numbers, value=1, command=lambda:self.update_display_flags("show_numbers", val=1))
     check_button_show_node_numbers.grid(row=2, column=0)

     # show linear node numbers: 
     check_button_show_node_linnumbers = Radiobutton(display_flags_frame, text="Node Indices (Linear)", variable=self.show_numbers, value=2, command=lambda:self.update_display_flags("show_numbers", val=2))
     check_button_show_node_linnumbers.grid(row=2, column=1)

     # show element numbers: 
     check_button_show_element_numbers = Radiobutton(display_flags_frame, text="Element Indices", variable=self.show_numbers, value=3, command=lambda:self.update_display_flags("show_numbers", val=3))
     check_button_show_element_numbers.grid(row=2, column=2)
     
     # show face numbers: 
     check_button_show_face_numbers = Radiobutton(display_flags_frame, text="Face Indices", variable=self.show_numbers, value=4, command=lambda:self.update_display_flags("show_numbers", val=4))
     check_button_show_face_numbers.grid(row=2, column=3)

     # show no numbers: 
     check_button_show_no_numbers = Radiobutton(display_flags_frame, text="No Indices", variable=self.show_numbers, value=0, command=lambda:self.update_display_flags("show_numbers", val=0))
     check_button_show_no_numbers.grid(row=2, column=4)
     check_button_show_no_numbers.select() # that has to match with the default value 1!

     # show springs: 
     check_button_show_springs = Checkbutton(display_flags_frame, text="Springs", variable=self.show_springs, command=lambda:self.update_display_flags("show_springs"))
     check_button_show_springs.grid(row=3, column=0)
     check_button_show_springs.select()

     # show pinned_nodes: 
     check_button_show_pinned = Checkbutton(display_flags_frame, text="Pinned Nodes", variable=self.show_pinned, command=lambda:self.update_display_flags("show_pinned"))
     check_button_show_pinned.grid(row=4, column=0)
     check_button_show_pinned.select() # that has to match with the default value 1!
 
     # Outer simulation box
     check_button_show_box = Radiobutton(display_flags_frame, text="Simulation Box (outline)", variable=self.show_box, value=1, command=lambda:self.update_display_flags("show_box", val=1))
     check_button_show_box.grid(row=5, column=0)
     check_button_show_box.select() # that has to match with the default value 1!     

     # Whole simulation box
     check_button_show_whole_box = Radiobutton(display_flags_frame, text="Simulation Box (whole)", variable=self.show_box, value=2, command=lambda:self.update_display_flags("show_box", val=2))
     check_button_show_whole_box.grid(row=5, column=1)

     # No simulation box
     check_button_show_box = Radiobutton(display_flags_frame, text="No box", variable=self.show_box, value=0, command=lambda:self.update_display_flags("show_box", val=0))
     check_button_show_box.grid(row=5, column=2)

     # load the trajectory:
     check_button_do_load_trajectory = Checkbutton(display_flags_frame, text="Load trajectory", variable=self.load_trajectory, command=lambda:self.update_display_flags("load_trajectory"))
     check_button_do_load_trajectory.grid(row=6, column=0)
     check_button_do_load_trajectory.select() # that has to match with the default value 1! 

     # flags
     self.animate = False
     self.display_window_exists = False
     self.there_is_something_to_send_to_display_window = False
     self.change_frame_to = -1
     	
     self.selected_index = 0
     self.selected_blob = 0
     self.selected_conformation = 0

     self.num_frames_to_read = float("inf")


    
 #################################################
  # # # # Update display_flags from buttons # # # 
 #################################################
  def update_display_flags(self, key, val=-1):

     # If unset (i.e. checkbutton)
     if val == -1:
	self.display_flags[key] = (self.display_flags[key] + 1) % 2
     else:
	self.display_flags[key] = val

     # WARNINGs:
     #NOT_IMPLEMENTED = ["show_element_numbers", "show_face_numbers"]
     #if NOT_IMPLEMENTED.count(key):
     #  print key, " functionality is still under development."



  # # # # # # # # # # # # # # # # # # # # # #
  # # Open dialogue for FFEA input file # # # 
  # # # # # # # # # # # # # # # # # # # # # # 
  def choose_ffea_file_to_load(self):
     # set up the options for the open file dialog box
     options = {}
     options['defaultextension'] = '.ffea'
     options['filetypes'] = [('ffea files', '.ffea'), ('all files', '.*')]
     options['initialdir'] = os.getcwd()
     options['title'] = 'Load ffea file'

     # Ask user to select a file
     ffea_fname = tkFileDialog.askopenfilename(**options)
     if len(ffea_fname) == 0:
             return

     # load the file
     self.load_ffea(ffea_fname)



  # # # # # # # # # # # # # # # # # # # # # #
  # # # # # # Load the FFEA file # # # # # # 
  # # # # # # # # # # # # # # # # # # # # # # 
  def load_ffea(self, ffea_fname):
  	
	# Try to reset previous system and update
	self.num_frames = 0
	self.num_loads += 1

	# Check if given file exists
	if os.path.isfile(ffea_fname) == False:
		print "No such file:", ffea_fname
		return
	else: 
		self.ffea_fname = ffea_fname
        
	print "Loading ffea file: " + self.ffea_fname
    
    # Load script (comments are now removed inside this module, by the way :) )
	self.script = FFEA_script.FFEA_script(self.ffea_fname)
	p = self.script.params
	bl = self.script.blob
	
	# See whether or not to remove traj file (make this better later i.e. rolling loading by storing file pointers)
	if self.display_flags['load_trajectory'] == 0:
		print "requested not to load the trajectory"
		p.trajectory_out_fname = None
        
    # Rebuild the script object depending on whether or not there is a trajectory (keep only first conformation)
	if p.trajectory_out_fname == None:
		for i in range(p.num_blobs):
			p.num_conformations[i] = 1
			bl[i].conformation = [bl[i].conformation[0]]     
    	
    #
    # Build the blob objects one at a time
    #
	self.blob_list = [[None for j in range(p.num_conformations[i])] for i in range(p.num_blobs)]
    
	idnum = 0
	for b in bl:
		bindex = bl.index(b)
		for c in b.conformation:
			cindex = b.conformation.index(c)
			ffea_id_string = "lol"
			print "\nLoading blob " + str(bindex) + ", conformation " + str(cindex)
			new_blob = Blob.Blob()
			#new_blob.load(blob_number, blob_index, conformation_index, blob_nodes[i], blob_top[i], blob_surface[i], blob_vdw[i], scale, blob_motion_state[i], blob_pin[i], blob_mat[i], blob_binding[i], blob_centroid_pos, blob_rotation, ffea_path)
			new_blob.load(idnum, bindex, cindex, self.script)
			new_blob.set_num_loads(self.num_loads)
     
			self.blob_list[bindex][cindex] = new_blob
			new_blob_name = ffea_id_string + "#" + str(bindex) + ", " + str(cindex)
			info_string = "Name:\t" + ffea_id_string + "\nConformation:\t" + str(cindex) + "\nNodes:\t" + c.nodes + "\nTopology:\t" + c.topology + "\nSurface:\t" + c.surface + "\nVdW:\t" + c.vdw + "\npin:\t" + c.pin + "\nMotion State:\t" + c.motion_state + "\n"
			add_blob_info = {'name': new_blob_name, 'info': info_string}
			
			idnum += 1
                 
    
	# Send binding sites to control
	binding_sites = [[0 for j in range(self.script.params.num_conformations[i])] for i in range(self.script.params.num_blobs)]
	for i in range(self.script.params.num_blobs):
		for j in range(self.script.params.num_conformations[i]):
			if self.blob_list[i][j].bsites != None:
				binding_sites[i][j] = self.blob_list[i][j].bsites.num_binding_sites

	# Rescale and translate initial system if necessary
	self.global_scale = 1e-10	# angstroms cos pymol works in angstroms and FFEA works in SI
	self.global_scale = 1.0 / self.global_scale

	# Rescale blobs
	for b in self.blob_list:
		for c in b:
			c.set_global_scale(self.global_scale)

	# Move simulation into box, if necessary
	world_centroid = np.array([0.0, 0.0, 0.0])
	shift = np.array([0.0, 0.0, 0.0])
	total_num_nodes = 0

	# Load all initial blobs and get a global centroid. Set secondary blobs to have placeholder 'None' frames
	for b in self.blob_list:

		for c in b:
			if b.index(c) == 0:
		
				c.set_nodes_as_frame()
				print c.frames[0].pos[0]
				x, y, z = c.get_centroid(0)
				world_centroid[0] += x * c.node.num_nodes
				world_centroid[1] += y * c.node.num_nodes
				world_centroid[2] += z * c.node.num_nodes
				total_num_nodes += c.node.num_nodes
			else:
				c.set_dead_frame()
			

	# Calculate global centroid
	world_centroid *= 1.0 / total_num_nodes	

	# Build box object if necessary
	if p.calc_vdw == 0:
		self.box = False
		self.box_x = 0.0
		self.box_y = 0.0
		self.box_z = 0.0
	else:
		try:

			# Do we need to calculate the box? Double the rounded up size of the system
			if p.es_N_x < 1 or p.es_N_y < 1 or p.es_N_z < 1:
				dims = self.get_system_dimensions(0)
				p.es_N_x = int(dims[0] * (p.kappa / (p.es_h * self.global_scale)))
				p.es_N_y = int(dims[1] * (p.kappa / (p.es_h * self.global_scale)))
				p.es_N_z = int(dims[2] * (p.kappa / (p.es_h * self.global_scale)))

			self.box = True
			self.box_x = 2 * (1.0 / p.kappa) * p.es_h * p.es_N_x
			self.box_y = 2 * (1.0 / p.kappa) * p.es_h * p.es_N_y
			self.box_z = 2 * (1.0 / p.kappa) * p.es_h * p.es_N_z
	
			# Does it exist? Realllllly?? If it's this small, it doesn't. OK?!!
			if self.box_x <= 1e-10 or self.box_y <= 1e-10 or self.box_z <= 1e-10:
				self.box = False
					
		except:
			self.box = False
			self.box_x = 0.0
			self.box_y = 0.0
			self.box_z = 0.0

		
	# Rescale box
	self.box_x *= self.global_scale
	self.box_y *= self.global_scale
	self.box_z *= self.global_scale

	shift[0] = self.box_x / 2.0 - world_centroid[0]
	shift[1] = self.box_y / 2.0 - world_centroid[1]
	shift[2] = self.box_z / 2.0 - world_centroid[2]

    
	# Shift all blobs if necessary
	for b in self.blob_list:
		if p.calc_vdw == 1 and p.move_into_box == 1:
			b[0].frames[0].translate(shift)
    		

    	# Now all blobs should have a single frame. Primary blobs should be in their starting configuration.
	# Secondary blobs should have a "None" placeholder. Therefore, we can draw it!
    		       
	# Now load trajectory (always run this function, regardless of stuff. It returns if anything is wrong)
	#if (p.trajectory_out_fname != None): # and (self.display_flags['load_trajectory'] == 1):
	self.load_trajectory_thread = threading.Thread(target=self.load_trajectory, args=(p.trajectory_out_fname, ))
	self.load_trajectory_thread.start()

	# Make sure we have at least 1 frame sorted before continuing, so main thread doesn't overtake
	#print self.num_frames, p.trajectory_out_fname
	#while(self.num_frames < 1):
	#	if p.trajectory_out_fname == None:
	#		# increase the frames to 1, so that the structure is displayed.
	#		self.num_frames = 1
	#		self.draw_stuff()
	#		break
	#	else:
	#		pass


  def load_trajectory(self, trajectory_out_fname):
	
	#
	# All blobs already have the first frame. They will keep this permanently.
	# All subsequent frames will be readed, loaded, drawn and deleted until failure
	#	

	# Load header and skip first frame (we already have it from the node files)
	traj = FFEA_trajectory.FFEA_trajectory(trajectory_out_fname, load_all = 0)
	try:
		failure = traj.skip_frame()
	except:
		failure = 1

	# Draw first frame
	self.num_frames = 1
	self.draw_frame(self.num_frames - 1)

	# If necessary, stop now (broken traj or user asked for)
	if traj.num_blobs == 0 or failure == 1 or self.display_flags['load_trajectory'] == 0:		
		return

	# Else, load rest of trajectory 1 frame at a time, drawing and deleting as we go
	while True:
		
		# Get frame from traj
		if traj.load_frame() == 0:
			
			# Scale traj frame
			traj.scale(self.global_scale, 0)

			# Load into blob objects asnd increment frame count
			self.add_frame_to_blobs(traj)
			self.num_frames += 1

			# Draw whole frame (if above worked, these should work no problem...)
			self.draw_frame(self.num_frames - 1)

			# Delete frames from memory
			traj.delete_frame()
			self.remove_frame_from_blobs()

		else:
			break
	return

  def get_system_dimensions(self, findex):
	maxdims = np.array([float("-inf"),float("-inf"),float("-inf")])	
	mindims = np.array([float("inf"),float("inf"),float("inf")])
	dims = np.array([0.0,0.0,0.0])

	try:
		for b in self.blob_list:
			for c in b:
				try:
					for p in c.frames[findex].pos:
						for i in range(3):
							if p[i] > maxdims[i]:
								maxdims[i] = p[i]
							if p[i] < mindims[i]:
								mindims[i] = p[i]

					
				except:
					continue
		
		for i in range(3):
			dims[i] = maxdims[i] - mindims[i]

		return dims
				
	except:
		return np.array([0.0,0.0,0.0])
	
  #def load_trajectory(self, trajectory_out_fname):
  #
  #	# This function will load the trajectory by:
  #		# Loading header.
  #		# Skip first frame (we already have it). 
  #		# Load frames 1 at a time and leave thread open to be manually activated by user and constantly check for newly written frames
  #		
  #		# Load header stuff automatically
  #		traj = FFEA_trajectory.FFEA_trajectory(trajectory_out_fname, load_all = 0)
  #		
  #		# Check for failure!
 # 		if traj.num_blobs == 0:
 # 		
 # 			# This will activate the draw_stuff for a single frame
 # 			print "Error. Problem with trajectory file. Cannot load."
 # 			self.script.params.trajectory_out_fname = None
#			return
#	
#		# Skip first frame as we already have it
#		traj.skip_frame()
#		
#		# Set num_frames for external stuff
#		self.num_frames = 1
#		
#		# Now, let's load a trajectory (while we can)
#		while True:
#		
#			# If user wants frames, give them frames
#			if self.display_flags['load_trajectory'] == 1:
#			
#				# Get a frame
#				if traj.load_frame() == 0:
#				
#					# Success! We got a new frame. Add it to blob
#					self.add_frame_to_blobs(traj)
#					self.num_frames += 1
#					self.draw_frame()
#					
#					# And clear the blob
#					traj.clear_frame()
#					self.remove_frame_from_blobs()
#				else:
#					
#					# All failures move to the beginning of what will be the next available frame. Wait a bit and continue
#					print self.num_frames
#					break
#					time.sleep(10)
#		
#			else:
#			
#				# Check again every 3 seconds
#				time.sleep(3)
#			
#			if self.num_frames > 1:
#				cmd.mset("1-"+str(self.num_frames))
#				if self.num_frames > 2:
#					cmd.mplay()
          
  def add_frame_to_blobs(self, traj, index = -1):
  	
  	for i in range(self.script.params.num_blobs):
  		for j in range(self.script.params.num_conformations[i]):
  			self.blob_list[i][j].frames.append(traj.blob[i][j].frame[index])
  			self.blob_list[i][j].num_frames += 1
  			
  def remove_frame_from_blobs(self, index = -1):
  	
  	for i in range(self.script.params.num_blobs):
  		for j in range(self.script.params.num_conformations[i]):
  			del self.blob_list[i][j].frames[index]
  			self.blob_list[i][j].num_frames -= 1
  		
  def init_vars(self):

	# num times loaded
	self.num_loads = 0

	# camera
	# self.orientation = Quaternion()
	self.z = 1

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

	self.display_flags = {
		'show_solid': 1, ## PYMOL OK
		'show_mesh': 0,
		'show_numbers': 0, ## PYMOL OK
		'show_pinned': 1,
		'show_shortest_edge': 0,
		'show_springs': 1,
		'show_box': 1,
		'load_trajectory': 1, ## PYMOL OK
		'show_inverted': 0,}

	self.buttons = {'show_mesh' : self.show_mesh,
		'show_solid' : self.show_solid,}

     # start the buttons to the default values
	for k in self.buttons.keys():
		self.buttons[k].set(self.display_flags[k])

	self.selected_index = 0
	self.selected_blob = 0
	self.selected_conformation = 0

	self.offset_x = 0
	self.offset_y = 0
	self.offset_z = 0

	# Assume box exists
	self.box = True
	self.box_x = -1
	self.box_y = -1
	self.box_z = -1
	self.springs = None

	self.modifying_frame = False

	self.recording = 0
	self.movie_dir = "__temp__FFEA_viewer_movie_dir__"

	self.projection = "perspective"


  def draw_frame(self, index):

	# Blobs should only ever have at most 2 frames on them, the initial one and the currently loaded one. So...
	frame_real_index = index

	if index > 0:
		frame_stored_index = 1
	else:
		frame_stored_index = 0
		
	# World first
	if self.display_flags['show_box'] != 0 and self.box == True:
		self.draw_box(frame_real_index)

	for i in range(self.script.params.num_blobs):
		for j in range(self.script.params.num_conformations[i]):
			self.blob_list[i][j].draw_frame(frame_stored_index, frame_real_index, self.display_flags)
  	
  def draw_stuff(self):

    # World first
    if self.display_flags['show_box'] == 1:
	self.draw_box()

    for f in range(self.num_frames):
      for i in range(self.script.params.num_blobs):
          for j in range(self.script.params.num_conformations[i]):
              self.blob_list[i][j].draw_frame(f, self.display_flags) ## PLUGIN OUT



      if self.springs != None and self.display_flags['show_springs'] == 1:
         self.draw_springs()

  def draw_box(self, f):
	
	# A cube has 8 vertices and 12 sides. A hypercube has 16 and 32! "Whoa, that's well cool Ben!" Yeah, ikr 
	obj = [BEGIN, LINES]
	
	# If only outline, no need to loop over entire plane
	if self.display_flags['show_box'] == 1:
		
		step_x = self.box_x
		step_y = self.box_y
		step_z = self.box_z

		# Loop over the three planes
		for i in range(2):
			for j in range(2):
				
				# Get a pair of vertices
				verts = [[i * step_x, j * step_y, 0.0], [i * step_x, j * step_y, self.box_z]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		for i in range(2):
			for j in range(2):
				
				# Get a pair of vertices
				verts = [[0.0, i * step_y, j * step_z], [self.box_x, i * step_y, j * step_z]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		for i in range(2):
			for j in range(2):
				
				# Get a pair of vertices
				verts = [[j * step_x, 0.0, i * step_z], [j * step_x, self.box_y, i * step_z]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])
			
	elif self.display_flags['show_box'] == 2:

		step_x = self.box_x / self.script.params.es_N_x
		step_y = self.box_y / self.script.params.es_N_y
		step_z = self.box_z / self.script.params.es_N_z

		# Loop over the three planes
		for i in range(self.script.params.es_N_x + 1):
			for j in range(self.script.params.es_N_y + 1):
				
				# Get a pair of vertices
				verts = [[i * step_x, j * step_y, 0.0], [i * step_x, j * step_y, self.box_z]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		for i in range(self.script.params.es_N_y + 1):
			for j in range(self.script.params.es_N_z + 1):
				
				# Get a pair of vertices
				verts = [[0.0, i * step_y, j * step_z], [self.box_x, i * step_y, j * step_z]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		for i in range(self.script.params.es_N_z + 1):
			for j in range(self.script.params.es_N_x + 1):
				
				# Get a pair of vertices
				verts = [[j * step_x, 0.0, i * step_z], [j * step_x, self.box_y, i * step_z]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

	#verts = [[0.0,0.0,0.0], [self.box_x,0.0,0.0], [self.box_x,0.0,self.box_z], [0.0,0.0,self.box_z], [0.0,self.box_y,0.0], [self.box_x,self.box_y,0.0], [self.box_x,self.box_y,self.box_z], [0.0,self.box_y,self.box_z]]
	
	#for i in range(4):
	#	obj.extend([VERTEX, verts[i][0], verts[i][1], verts[i][2]])
	#	obj.extend([VERTEX, verts[(i + 1) % 4][0], verts[(i + 1) % 4][1], verts[(i + 1) % 4][2]])
#
#		obj.extend([VERTEX, verts[i][0], verts[i][1], verts[i][2]])
#		obj.extend([VERTEX, verts[i + 4][0], verts[i + 4][1], verts[i + 4][2]])
#
#		obj.extend([VERTEX, verts[i + 4][0], verts[i + 4][1], verts[i + 4][2]])
#		obj.extend([VERTEX, verts[(i + 1) % 4 + 4][0], verts[(i + 1) % 4 + 4][1], verts[(i + 1) % 4 + 4][2]])

				
					

	obj.append(END)
	cmd.load_cgo(obj, "Simulation Box", f)

  def draw_springs(self):

      for s in self.springs.spring:

         # Get correct frames
         correct_frame = [-1 for i in range(self.script.params.num_blobs)]
         for i in range(self.script.params.num_blobs):
            if self.blob_list[i][0].state == "STATIC":
               correct_frame[i] = 0
         # print "correct_frame: ", correct_frame

         # Draw, because this spring exists
         springjoints = np.array([self.blob_list[s.blob_index[i]][s.conformation_index[i]].frames[correct_frame[s.blob_index[i]]].node_list[s.node_index[i]][0:3] for i in range(2)])

         # Axes for helix
         zax = springjoints[1] - springjoints[0]
         xax = np.cross(zax, np.array([1.0,0]))
         yax = np.cross(zax, xax)

         xax = xax / np.linalg.norm(xax)
         yax = yax / np.linalg.norm(yax)

         l = np.linalg.norm(zax)

         zax = zax / l

         # Radius of helix (let original radius be 5A, poisson ration = 0.01)
         r = 5 - 0.01 * (l - s.l)

         # We want 5 spins, say, so pitch:
         c = l / (10 * np.pi)

         # Draw 40 nodes. Equation is r = r0 + (Rcos(t), Rsin(t), ct)
         step = (10 * np.pi) / 40
         obj = [ BEGIN, LINES, LINEWIDTH, 4.0 ]
         # obj.extend( [ COLOR, 192/255.0, 192/255.0, 192/255.0 ] )
         for i in range(40):
            tstart = step * i
            tend = step * (i + 1)
            verts = springjoints[0] + np.array([r * np.cos(tstart) * xax[i] + r * np.sin(tstart) * yax[i] + c * tstart * zax[i] for i in range(3)])
            obj.extend( [ VERTEX, verts[0], verts[1], verts[2] ] )
            verts = springjoints[0] + np.array([r * np.cos(tend) * xax[i] + r * np.sin(tend) * yax[i] + c * tend * zax[i] for i in range(3)])
            obj.extend( [ VERTEX, verts[0], verts[1], verts[2] ] )

         obj.append(END)
         cmd.load_cgo(obj, "string_" + str(self.springs.spring.index(s)), max(correct_frame))
         


