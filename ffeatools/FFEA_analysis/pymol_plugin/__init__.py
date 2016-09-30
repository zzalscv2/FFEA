import sys, os, time
import numpy as np

from pymol import cmd
from pymol.callback import Callback
import warnings

import warnings

try:
    from mtTkinter import *
except ImportError:
    warnings.warn("HORRIBLE DANGER: Tkinter is not thread-safe. Viewer will now crash. Please install mtTKinter.", RuntimeWarning)
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
import FFEA_turbotrajectory
import FFEA_surface

from numpy.random import randint as rint

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

     self.root.geometry("680x200")

     self.root.title("FFEA")

     top_frame = Frame(self.root)
     top_frame.pack()

     menubar = Menu(top_frame)

     filemenu = Menu(menubar, tearoff=0)
     filemenu.add_command(label="Load 'ffea' file", command=self.choose_ffea_file_to_load)
     menubar.add_cascade(label="File", menu=filemenu)
     self.root.config(menu=menubar)

     # PLUGIN (separated into mutually exclusive sets. Devs take note!)



     self.init_vars()
     self.system_name = StringVar(self.root, value=self.display_flags['system_name'])
     self.do_load_trajectory = IntVar(self.root, value=self.display_flags['load_trajectory'])
     self.show_box = IntVar(self.root, value=self.display_flags['show_box'])
     self.show_numbers = IntVar(self.root, value=self.display_flags['show_numbers'])
     self.show_pinned = IntVar(self.root, value=self.display_flags['show_pinned'])
     self.show_vdw = IntVar(self.root, value=self.display_flags['show_vdw'])
     self.show_springs = IntVar(self.root, value=self.display_flags['show_springs'])
     self.show_solid = IntVar(self.root, value=self.display_flags['show_solid'])
     self.matparam = StringVar(self.root, value=self.display_flags['matparam'])
     self.show_mesh = IntVar(self.root, value=self.display_flags['show_mesh'])
     self.show_shortest_edge = IntVar(self.root, value=self.display_flags['show_shortest_edge'])


 
     # # Display flags frame
     display_flags_frame = Frame(self.root)
     display_flags_frame.pack(anchor=CENTER, expand=True)


     # propose a system name:
     label_system_name = Label(display_flags_frame, text="System name:")
     label_system_name.grid(row=0, column=0, sticky=E)
     text_button_system_name = Entry(display_flags_frame, text="load as:", textvariable=self.system_name, validate="focus", validatecommand=lambda:self.update_display_flags("system_name", val=-2, text=self.system_name.get()))
     text_button_system_name.grid(row=0, column=1, sticky=W)
     
     random_name_button = Button(display_flags_frame, text="Random Name", command=lambda:self.new_system_name());
     random_name_button.grid(row=0, column=2, sticky=W)

     label_display = Label(display_flags_frame, text="Display:")
     label_display.grid(row=1, column=0, sticky=E)

     # show springs: 
     check_button_show_springs = Checkbutton(display_flags_frame, text="Springs", variable=self.show_springs, command=lambda:self.update_display_flags("show_springs"))
     check_button_show_springs.grid(row=1, column=1, sticky=W)


     # show pinned_nodes: 
     check_button_show_pinned = Checkbutton(display_flags_frame, text="Pinned Nodes", variable=self.show_pinned, command=lambda:self.update_display_flags("show_pinned"))
     check_button_show_pinned.grid(row=1, column=2, sticky=W)
 
    # # show vdw
     check_button_show_vdw = Checkbutton(display_flags_frame, text="VdW Faces", variable=self.show_vdw, command=lambda:self.update_display_flags("show_vdw"))
     check_button_show_vdw.grid(row=1, column=3, sticky=W)

     # # show solid:

     label_solid = Label(display_flags_frame, text="Show Solid:")
     label_solid.grid(row=2, column=0, sticky=E)

     SolidModes = [(0, "Plain Solid", 1),\
                   (1, "Material", 2),\
                   (3, "No Solid", 0)] 
     for col, text, mode in SolidModes:
        check_button_show_solid = Radiobutton(display_flags_frame, text=text, variable=self.show_solid, value=mode, command=lambda:self.update_display_flags("show_solid", val=self.show_solid.get()))
        check_button_show_solid.grid(row=2, column=col+1, sticky=W) # first col is label

     # Selectable box for material param
     spinbox_material_param = Spinbox(display_flags_frame, textvariable=self.matparam, values=("Density", "Shear Viscosity", "Bulk Viscosity", "Shear Modulus", "Bulk Modulus"), validate="focus", validatecommand=lambda:self.update_display_flags("matparam", val=-2, text=self.matparam.get()))
     spinbox_material_param.grid(row=2, column=3, sticky=W)

     # # show mesh:

     label_mesh = Label(display_flags_frame, text="Show Mesh:")
     label_mesh.grid(row=3, column=0, sticky=E)

     MeshModes = [(0, "Surface Mesh", 2),\
                  (1, "Whole Mesh", 1),\
                  (2, "No Mesh", 0)]
     for col, text, mode in MeshModes:
        check_button_show_mesh = Radiobutton(display_flags_frame, text=text, variable=self.show_mesh, value=mode, command=lambda:self.update_display_flags("show_mesh", val=self.show_mesh.get()))
        check_button_show_mesh.grid(row=3, column=col+1, sticky=W) # first col is label

     label_mesh= Label(display_flags_frame, text="Show Indices:")
     label_mesh.grid(row=4, column=0, sticky=E)
     
     # # show Numbers:

     self.show_numbers = StringVar(display_flags_frame)
     self.show_numbers.set("No Indices")
     index_option = OptionMenu(display_flags_frame, self.show_numbers, "Node Indices", "Node Indices (Linear)", "Element Indicies", "Face Indices", "No Indices", command=lambda x:self.update_display_flags("show_numbers", val=self.show_numbers.get()) )
     index_option.grid(row=4, column=1, sticky=W)
     
     button_node_indices_psuedoatoms = Button(display_flags_frame, text="Add node psuedoatoms", command=lambda:self.add_node_psuedoatoms())
     button_node_indices_psuedoatoms.grid(row=4, column=2, sticky=W)

     label_box= Label(display_flags_frame, text="Show Box:")
     label_box.grid(row=5, column=0, sticky=E)
 
     # Outer simulation box
     BoxModes = [(0, "Simulation Box (outline)", 1),\
                 (1, "Simulation Box (whole)", 2),\
                 (2, "No Box", 0)]
     for col, text, mode in BoxModes:
       check_button_show_box = Radiobutton(display_flags_frame, text=text, variable=self.show_box, value=mode, command=lambda:self.update_display_flags("show_box", val=self.show_box.get()))
       check_button_show_box.grid(row=5, column=col+1, sticky=W) # first col is label

     label_traj= Label(display_flags_frame, text="Load:")
     label_traj.grid(row=6, column=0, sticky=E)

     ## # Trajectory Radiobutton # #
     TrjModes = [("Trajectory", 1), \
                 ("System (Into box)", 2), \
                 ("System (Plainly)",3), \
                 ("CGO", 4)]
     for text, mode in TrjModes:
       check_button_do_load_trajectory = Radiobutton(display_flags_frame, text=text, variable=self.do_load_trajectory, value=mode, command=lambda:self.update_display_flags("load_trajectory", val=self.do_load_trajectory.get()))
       check_button_do_load_trajectory.grid(row=6, column=mode, sticky=W)
  
     label_traj= Label(display_flags_frame, text="Highlight element(s):")
     label_traj.grid(row=7, column=0, sticky=E)
     
     self.highlight = StringVar(self.root, value=self.display_flags['highlight'])
     highlight_entry_box = Entry(display_flags_frame, textvariable=self.highlight, validate="focus", validatecommand=lambda:self.update_display_flags("highlight", val=-2, text=self.highlight.get()))
     highlight_entry_box.grid(row=7, column=1, sticky=W)

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
  # # use val = -2 for strings (Entries)
  # #     val = -1 for binary choices (Checkboxes)
  # #     val > 0, integer, for Radiobuttons 
 #################################################
  def update_display_flags(self, key, val=-1, text=""):

     # If unset (i.e. checkbutton)
     if val == -2:
       self.display_flags[key] = text
       return True
     elif val == -1:
       self.display_flags[key] = (self.display_flags[key] + 1) % 2
     else:
       self.display_flags[key] = val

     print key, ": ", self.display_flags[key]


  def new_system_name(self):

	self.system_name.set(self.system_names[rint(0, len(self.system_names) - 1)])
	self.display_flags["system_name"] = self.system_name.get()

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
#     self.root.destroy()
     self.load_ffea(ffea_fname)



  # # # # # # # # # # # # # # # # # # # # # #
  # # # # # # Load the FFEA file # # # # # # 
  # # # # # # # # # # # # # # # # # # # # # # 
  def load_ffea(self, ffea_fname):
  	
	# Update display flags patch (the .get() function got the old spinbox value, so here it's definitely updated)
	self.display_flags['matparam'] = self.matparam.get()

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
	if (self.script.params == None):
		print "Something went wrong reading the FFEA input file", self.ffea_fname
		return
	p = self.script.params
	bl = self.script.blob
	
	# See whether or not to remove traj file (make this better later i.e. rolling loading by storing file pointers)
	if self.display_flags['load_trajectory'] == 2 or 3:
		print "Requested not to load the trajectory"
		#p.trajectory_out_fname = None
	if self.display_flags['load_trajectory'] == 3:
		print "Requested to show the coordinates as they are in the .node(s) file(s)"
		print "... equivalently, setting < move_into_box = 0 >"
		p.move_into_box = 0
        
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
	bindex = -1
	for b in bl:
		bindex += 1
		cindex = -1
		for c in b.conformation:
			cindex += 1
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
                 
    
	# Load some springs
	if self.display_flags['show_springs'] == 1:
		self.springs = FFEA_springs.FFEA_springs(self.script.spring)

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
	bindex = -1	
	for b in self.blob_list:
		bindex += 1
		cindex = -1
		for c in b:
			cindex += 1
			if cindex == 0:
		
				c.set_nodes_as_frame()

				x, y, z = c.calc_centroid(0)
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
		self.box_exists = False
		self.box = np.array([0.0,0.0,0.0])
	else:
		try:

			# Do we need to calculate the box? Double the rounded up size of the system
			for i in range(3):
				if p.es_N[i] < 1:
					dims = self.get_system_dimensions(0)
					for j in range(3):
						p.es_N[j] = 2 * int(np.ceil(dims[j] * (p.kappa / (p.es_h * self.global_scale))))
					break

			self.box_exists = True
			self.box = (1.0 / p.kappa) * p.es_h * p.es_N
	
			# Does it exist? Realllllly?? If it's this small, it doesn't. OK?!!
			for i in self.box:
				if i <= 1e-10:
					self.box_exists = False
					break

		except:
			self.box_exists = False
			self.box = np.array([0.0,0.0,0.0])

		
	# Rescale box
	self.box *= self.global_scale

	# Shift all blobs to center of box if necessary
	shift = 0.5 * self.box - world_centroid
	if p.calc_vdw == 1 and p.move_into_box == 1:
		for b in self.blob_list:
			b[0].frames[0].translate(shift)
    		

	# Now, apply PBC if necessary
	if p.calc_vdw == 1:
		for b in self.blob_list:
			trans = np.array([0.0,0.0,0.0])
			cent = b[0].frames[0].calc_centroid()
			print "Centroid = ", cent
			if self.box_exists:
				for i in range(3):
					if cent[i] > self.box[i]:
						trans[i] = -1 * self.box[i]
					elif cent[i] < 0:
						trans[i] = self.box[i]

				b[0].frames[0].translate(trans)
				print "Translation = ", trans

    	# Now all blobs should have a single frame. Primary blobs should be in their starting configuration.
	# Secondary blobs should have a "None" placeholder. Therefore, we can draw it!
    		       
	# Now load trajectory (always run this function, regardless of stuff. It returns if anything is wrong)
	#if (p.trajectory_out_fname != None): # and (self.display_flags['load_trajectory'] == 1):
	traj_fname = self.script.params.trajectory_out_fname
	cgo_fname = traj_fname.split(".")[0]+"_cgo.npy"
	cgo_index_fname = traj_fname.split(".")[0]+"_cgoindex.npy"
	if self.display_flags['load_trajectory'] == 4:
		if os.path.isfile(cgo_fname) == False:
			print("No cached traj found at "+cgo_fname+", generating one...")
			turbotraj = FFEA_turbotrajectory.FFEA_turbotrajectory()
			turbotraj.populate_turbotraj_from_ftj(self.script.params.trajectory_out_fname)
			turbotraj.create_cgo(self.script, self.display_flags)
			turbotraj.dump_cgo()
		self.load_cgo(cgo_fname, cgo_index_fname)
		#cmd.load_cgo(turbotraj.cgo, self.display_flags['system_name'], frame)
	else:
		self.load_trajectory_thread = threading.Thread(target=self.load_trajectory, args=(p.trajectory_out_fname, ))
		self.load_trajectory_thread.start()
		#self.load_trajectory(p.trajectory_out_fname)

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

	# Choose new system name for next run
	self.system_index += 1
	self.system_name = self.system_names[self.system_index]

  def get_normal(self, node0, node1, node2):
	ax = node1[0] - node0[0]
	ay = node1[1] - node0[1]
	az = node1[2] - node0[2]
	bx = node2[0] - node1[0]
	by = node2[1] - node1[1]
	bz = node2[2] - node1[2]

	return [az * by - ay * bz, ax * bz - az * bx, ay * bx - ax * by]

  def load_cgo(self, cgo_fname, cgo_index_fname):
      cgo = np.load(cgo_fname)
      cgo_index = np.load(cgo_index_fname)
      print("Loading the cgo object...")
      for frame in range(len(cgo_index)):
          cmd.load_cgo(cgo[frame], cgo_index[frame][0], str(cgo_index[frame][1]))

  def load_turbotrajectory(self, turbotraj):
      
    def setup(self, turbotraj):
        frames = range(len(turbotraj.turbotraj[0][0]))
        surfs = []

        # cerate a list of surfaces, one for each blob
        for i in range(len(self.blob_list)):
            surfs.append(self.script.load_surface(i)) 
        return surfs, frames
        
    def get_nodes_in_face(turbotraj, face):
        return [turbotraj.turbotraj[blob_num][0][frame][face.n[0]], turbotraj.turbotraj[blob_num][0][frame][face.n[1]], turbotraj.turbotraj[blob_num][0][frame][face.n[2]]]
    
    surfs, frames = setup(self, turbotraj)

    # for every frame, create a cgo object
    for frame in frames:
        print("Loading frame "+str(frame)+"...")
        sol = [ BEGIN, TRIANGLES ]

        # for each face in each surf, load the nodes into the cgo as triangles
        for blob_num in range(len(surfs)):
            for face in surfs[blob_num].face:
                nodexyz = get_nodes_in_face(turbotraj, face)
                norm = self.get_normal(nodexyz[0], nodexyz[1], nodexyz[2])
                sol.extend( [ NORMAL, -norm[0], -norm[1], -norm[2], VERTEX, nodexyz[0][0]*1000000000, nodexyz[0][1]*1000000000, nodexyz[0][2]*1000000000, VERTEX, nodexyz[1][0]*1000000000, nodexyz[1][1]*1000000000, nodexyz[1][2]*1000000000, VERTEX, nodexyz[2][0]*1000000000, nodexyz[2][1]*1000000000, nodexyz[2][2]*1000000000 ] )
        sol.append(END)#
        cmd.load_cgo(sol, self.display_flags['system_name'], frame)

    
    # Each trajectory is composed of several blobs, do it for all blobs
    # if each blob has several conformations, do it for all of those
    # skip when 'none' obviously
    # in each blob->conformation, consult the surface file
    # for each surf.face.n, grab the points at that index and draw a trinagle with them
    return
    
  def add_node_psuedoatoms(self):
      node_object_list = []
      for blob_num in range(len(self.blob_list)):
          traj = self.script.load_trajectory(1)
          node_object_list.append(traj.blob[blob_num])
      for node_object in range(len(node_object_list)):
          for conformation in node_object_list[node_object]:
              for node in range(len(conformation.frame[0].pos)):
                  if conformation.frame[0].pos[node] !=None:
                      cmd.pseudoatom(pos = (conformation.frame[0].pos[node]*10000000000).tolist(), name = str(node), color="black")
                          
  def add_node_psuedoatoms_from_nodes(self):
      node_object_list = []
      for blob_num in range(len(self.blob_list)):
          node_object_list.append(self.script.load_node(blob_num))
      for node_object in range(len(node_object_list)):
          for node in range(len(node_object_list[node_object].pos)):
              cmd.pseudoatom(pos = node_object_list[node_object].pos[node].tolist())   
        
 
  def load_trajectory(self, trajectory_out_fname):
	
	#
	# All blobs already have the first frame. They will keep this permanently.
	# All subsequent frames will be readed, loaded, drawn and deleted until failure
	#	

	# Load header and skip first frame (we already have it from the node files)
	try:
		traj = FFEA_trajectory.FFEA_trajectory(trajectory_out_fname, load_all = 0)
		try:
			failure = traj.skip_frame()
		except:
			failure = 1
	except(IOError):
		failure = 1	

	# Draw first frame
	self.num_frames = 1
	self.draw_frame(self.num_frames - 1)

	# If necessary, stop now (broken traj or user asked for)
	if failure == 1 or self.display_flags['load_trajectory'] != 1 or traj.num_blobs == 0:		
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

	# Finally show the "progress bar":
	if self.num_frames > 1:
		cmd.mset("1-"+str(self.num_frames))
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

	self.system_index = 0
	self.system_names = []
	
	# Change to any file of names you like
	fname = os.path.dirname(os.path.realpath(__file__)) + "/system_names_dbzcharacters.txt"
	with open(fname, "r") as f:
		for line in f:
			self.system_names.append(line.strip())

	self.display_flags = {
		'show_solid': 1, ## PYMOL OK
		'matparam': "Density",
		'show_mesh': 0,
		'show_numbers': 0, ## PYMOL OK
		'show_pinned': 1,
		'show_vdw': 0,
		'show_shortest_edge': 0,
		'show_springs': 1,
		'show_box': 0,
		'load_trajectory': 1, ## PYMOL OK
		'show_inverted': 0,
		'highlight': '',
      'system_name': self.system_names[rint(0, len(self.system_names) - 1)]}

	self.selected_index = 0
	self.selected_blob = 0
	self.selected_conformation = 0

	self.offset_x = 0
	self.offset_y = 0
	self.offset_z = 0

	# Assume box exists
	self.box_exists = True
	self.box = np.array([-1.0,-1.0,-1.0])
	self.springs = None

	self.modifying_frame = False

	self.recording = 0
	self.movie_dir = "__temp__FFEA_viewer_movie_dir__"

	self.projection = "perspective"


  def get_system_centroid(self, frameIndex = -1):

	cent = np.array([0.0,0.0,0.0])
	total_num_nodes = 0
	bindex = -1
	for b in self.blob_list:
		bindex += 1
		cindex = -1		
		for c in b:
			cindex += 1
			if cindex == 0:

				x, y, z = c.calc_centroid(frameIndex)
				cent[0] += x * c.node.num_nodes
				cent[1] += y * c.node.num_nodes
				cent[2] += z * c.node.num_nodes
				total_num_nodes += c.node.num_nodes
	cent *= 1.0 / total_num_nodes
	return cent

  def draw_frame(self, index):

	# Blobs should only ever have at most 2 frames on them, the initial one and the currently loaded one. So...
	frame_real_index = index

	if index > 0:
		frame_stored_index = 1
	else:
		frame_stored_index = 0
		
	# World first
	if self.display_flags['show_box'] != 0 and self.box_exists == True:
		self.draw_box(frame_real_index)

	if self.display_flags['show_springs'] == 1:
		self.draw_springs(frame_real_index)

	for i in range(self.script.params.num_blobs):
		for j in range(self.script.params.num_conformations[i]):
			self.blob_list[i][j].draw_frame(frame_stored_index, frame_real_index, self.display_flags)

  def draw_box(self, f):
	
	# A cube has 8 vertices and 12 sides. A hypercube has 16 and 32! "Whoa, that's well cool Ben!" Yeah, ikr 
	obj = [BEGIN, LINES]
	
	# If only outline, no need to loop over entire plane
	if self.display_flags['show_box'] == 1:
		
		step = self.box

		# Loop over the three planes
		for i in range(2):
			for j in range(2):
				
				# Get a pair of vertices
				verts = [[i * step[0], j * step[1], 0.0], [i * step[0], j * step[1], self.box[2]]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		for i in range(2):
			for j in range(2):
				
				# Get a pair of vertices
				verts = [[0.0, i * step[1], j * step[2]], [self.box[0], i * step[1], j * step[2]]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		for i in range(2):
			for j in range(2):
				
				# Get a pair of vertices
				verts = [[j * step[0], 0.0, i * step[2]], [j * step[0], self.box[1], i * step[2]]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])
			
	elif self.display_flags['show_box'] == 2:

		for i in range(3):
			step[i] = self.box[i] / self.script.params.es_N[i]

		# Loop over the three planes
		for i in range(self.script.params.es_N[0] + 1):
			for j in range(self.script.params.es_N[1] + 1):
				
				# Get a pair of vertices
				verts = [[i * step[0], j * step[1], 0.0], [i * step[0], j * step[1], self.box[2]]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		for i in range(self.script.params.es_N[1] + 1):
			for j in range(self.script.params.es_N[2] + 1):
				
				# Get a pair of vertices
				verts = [[0.0, i * step[1], j * step[2]], [self.box[0], i * step[1], j * step[2]]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		for i in range(self.script.params.es_N[2] + 1):
			for j in range(self.script.params.es_N[0] + 1):
				
				# Get a pair of vertices
				verts = [[j * step[0], 0.0, i * step[2]], [j * step[0], self.box[1], i * step[2]]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

				
					

	obj.append(END)
	cmd.load_cgo(obj, self.display_flags['system_name'] +"_Simulation_Box", f + 1)

  def draw_springs(self, f):

      for s in self.springs.spring:

         # Get correct frames
         correct_frame = [-1 for i in range(self.script.params.num_blobs)]
         for i in range(self.script.params.num_blobs):
            if self.blob_list[i][0].motion_state == "STATIC":
               correct_frame[i] = 0
         # print "correct_frame: ", correct_frame

         # Draw, because this spring exists
         springjoints = np.array([self.blob_list[s.blob_index[i]][s.conformation_index[i]].frames[correct_frame[s.blob_index[i]]].pos[s.node_index[i]][0:3] for i in range(2)])

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
         cmd.load_cgo(obj, self.display_flags['system_name'] + "_string_" + str(self.springs.spring.index(s)), f + 1)
         


