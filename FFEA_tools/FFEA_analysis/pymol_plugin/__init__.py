import sys, os
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

# FFEA stuff
import FFEA_script
import FFEA_trajectory

# do Ben's springs:
import FFEA_springs

# PyMOL stuff:
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain



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
     self.root.geometry("155x170")
     self.root.title("FFEA")

     top_frame = Frame(self.root)
     top_frame.pack()

     menubar = Menu(top_frame)

     filemenu = Menu(menubar, tearoff=0)
     filemenu.add_command(label="Load 'ffea' file", command=self.choose_ffea_file_to_load)
     menubar.add_cascade(label="File", menu=filemenu)
     self.root.config(menu=menubar)

     # PLUGIN
     self.show_mesh = IntVar()
     self.show_mesh_surf = IntVar()
     self.show_solid = IntVar()
     self.show_node_numbers = IntVar()
     self.show_element_numbers = IntVar()
     self.show_face_numbers = IntVar()
     self.show_springs = IntVar()
     self.do_load_trajectory = IntVar()
     self.init_vars()
 
     # # Display flags frame
     display_flags_frame = Frame(self.root, relief=SUNKEN, bd=1)
     # display_flags_frame.place(relx=.5, rely=.75, anchor="c")
     display_flags_frame.pack(anchor=CENTER, expand=True)

     # show mesh:
     check_button_show_mesh = Checkbutton(display_flags_frame, text="Mesh", variable=self.show_mesh, command=lambda:self.update_display_flags("show_mesh"))
     check_button_show_mesh.pack(anchor=W)

     # show mesh surf:
     check_button_show_mesh_surf = Checkbutton(display_flags_frame, text="Mesh surf", variable=self.show_mesh_surf, command=lambda:self.update_display_flags("show_mesh_surf"))
     check_button_show_mesh_surf.pack(anchor=W)

     # show solid: 
     check_button_show_solid = Checkbutton(display_flags_frame, text="Solid", variable=self.show_solid, command=lambda:self.update_display_flags("show_solid"))
     check_button_show_solid.pack(anchor=W)
     check_button_show_solid.select() # that has to match with the default value 1! 

     # show springs: 
     check_button_show_springs = Checkbutton(display_flags_frame, text="Springs", variable=self.show_springs, command=lambda:self.update_display_flags("show_springs"))
     check_button_show_springs.pack(anchor=W)
     
     # show node numbers: 
     check_button_show_node_numbers = Checkbutton(display_flags_frame, text="Node numbers", variable=self.show_node_numbers, command=lambda:self.update_display_flags("show_node_numbers"))
     check_button_show_node_numbers.pack(anchor=W)

     # show element numbers: 
     check_button_show_element_numbers = Checkbutton(display_flags_frame, text="Element numbers", variable=self.show_element_numbers, command=lambda:self.update_display_flags("show_element_numbers"))
     check_button_show_element_numbers.pack(anchor=W)
     
     # show face numbers: 
     check_button_show_face_numbers = Checkbutton(display_flags_frame, text="Face numbers", variable=self.show_face_numbers, command=lambda:self.update_display_flags("show_face_numbers"))
     check_button_show_face_numbers.pack(anchor=W)
     
     # load the trajectory:
     check_button_do_load_trajectory = Checkbutton(display_flags_frame, text="Load trajectory", variable=self.load_trajectory, command=lambda:self.update_display_flags("load_trajectory"))
     check_button_do_load_trajectory.pack(anchor=W)
     check_button_do_load_trajectory.select() # that has to match with the default value 1! 

     # flags
     self.animate = False
     self.display_window_exists = False
     self.there_is_something_to_send_to_display_window = False
     self.change_frame_to = -1
     	
     self.num_blobs = 0
     self.num_conformations = []
     self.selected_index = 0
     self.selected_blob = 0
     self.selected_conformation = 0

     self.num_frames_to_read = float("inf")


    
 #################################################
  # # # # Update display_flags from buttons # # # 
 #################################################
  def update_display_flags(self, key):

     # update display flags for key
     if self.display_flags[key] == 0:
       self.display_flags[key] = 1
     else:
       self.display_flags[key] = 0
     # print key, self.display_flags[key] 


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
        	bl[i] = [bl[i][0]]     
    
    # Build box object
    self.box_x = (1.0 / p.kappa) * p.es_h * p.es_N_x
    self.box_y = (1.0 / p.kappa) * p.es_h * p.es_N_y
    self.box_z = (1.0 / p.kappa) * p.es_h * p.es_N_z
    
    #
    # Build the blob objects one at a time
    #
    self.blob_list = [[None for j in range(p.num_conforamtions[i])] for i in range(p.num_blobs)]
    
    idnum = 0
   	for b in bl:
   	    bindex = bl.index(b)
   	    for c in b:
   	        cindex = b.index(c)

            print "\nLoading blob " + str(bindex) + ", conformation " + str(cindex)
            new_blob = Blob.Blob()
            #new_blob.load(blob_number, blob_index, conformation_index, blob_nodes[i], blob_top[i], blob_surface[i], blob_vdw[i], scale, blob_motion_state[i], blob_pin[i], blob_mat[i], blob_binding[i], blob_centroid_pos, blob_rotation, ffea_path)
            new_blob.load(idnum, bindex, cindex, self.script)     
                 self.blob_list[bindex][cindex] = new_blob
                 new_blob_name = ffea_id_string + "#" + str(bindex) + ", " + str(cindex)
                 info_string = "Name:\t" + ffea_id_string + "\nConformation:\t" + str(cindex) + "\nNodes:\t" + c.nodes + "\nTopology:\t" + c.topology + "\nSurface:\t" + c.surface + "\nVdW:\t" + c.vdw + "\npin:\t" + c.pin + "\nMotion State:\t" + b.motion_state + "\n"
                 add_blob_info = {'name': new_blob_name, 'info': info_string}
                 # self.speak_to_control.send({'add_blob': add_blob_info}) ## PLUGIN OUT 
                 idnum += 1
                 
    
    # Rescale and translate initial system if necessary
    # Send binding sites to control
    binding_sites = [[0 for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]
    for i in range(self.num_blobs):
        for j in range(self.num_conformations[i]):
        	if self.blob_list[i][j].num_binding_sites != None:
                binding_sites[i][j] = self.blob_list[i][j].bsites.num_binding_sites

     # Rescale and translate initial system if necessary
     global_scale = float("inf")
     for blob in self.blob_list:
         if blob[0].scale < global_scale:
             global_scale = blob[0].scale

     global_scale = 1.0 / global_scale

     # Rescale box
     self.box_x *= global_scale
     self.box_y *= global_scale
     self.box_z *= global_scale

     # Rescale blobs
     for b in self.blob_list:
         for c in b:
             c.set_global_scale(global_scale)

     # Move simulation into box, if necessary
     world_centroid = np.array([0.0, 0.0, 0.0])
     shift = np.array([0.0, 0.0, 0.0])
     total_num_nodes = 0

     # Load STATIC blobs and get a global centroid
     for b in self.blob_list:

         b[0].set_nodes_as_frame()
         x, y, z = b[0].get_centroid(0)
         world_centroid[0] += x * b[0].node.num_nodes
         world_centroid[1] += y * b[0].node.num_nodes
         world_centroid[2] += z * b[0].node.num_nodes
         total_num_nodes += b[0].node.num_nodes

     world_centroid *= 1.0 / total_num_nodes	
     	
     shift[0] = self.box_x / 2.0 - world_centroid[0]
     shift[1] = self.box_y / 2.0 - world_centroid[1]
     shift[2] = self.box_z / 2.0 - world_centroid[2]
     		
     # Shift all blobs if STATIC, or if there is no trajectory; clear frame if not
     for b in self.blob_list:
         if trajectory_out_fname == None:
             if self.calc_vdw == 1 and self.move_into_box == 1:
                 # b[0].frames[0].translate(shift)
                 print "--- not translated by: ", shift
         elif blob[0].state == "STATIC":
             if self.calc_vdw == 1 and self.move_into_box == 1:
                 b[0].frames[0].translate(shift)
                 print "--- translated by: ", shift
         else:
             b[0].frames = []
             b[0].num_frames = 0
             
     # Now load trajectory
     if (trajectory_out_fname != None): # and (self.display_flags['load_trajectory'] == 1):
         self.load_trajectory_thread = threading.Thread(target=self.load_trajectory, args=(trajectory_out_fname, ))
         self.load_trajectory_thread.start()

     ## Hold on calculating dimensions until at least one frame has been calculated from a trajectory, if it exists
     while(self.num_frames < 1):
       if trajectory_out_fname == None:
            # increase the frames to 1, so that the structure is displayed.
            self.num_frames = 1
            self.draw_stuff()
            break
       else:
           pass

     # Reset initial camera (dependent upon structure size)
     # dims = self.get_system_dimensions()
     # self.dimensions = [dims[i][1] - dims[i][0] for i in range(3)]
# 
#      if (self.dimensions[2] > self.dimensions[1]) and (self.dimensions[2] > self.dimensions[0]):
#          self.z = 2 * self.dimensions[2]
#      elif self.dimensions[0] > self.dimensions[1]:
#          self.z = self.dimensions[0] / (2 * np.tan(np.pi / 6.0))
#      else:
#         self.z = self.dimensions[1] / (2 * np.tan(np.pi / 6.0))

  def load_trajectory(self, trajectory_out_fname):
	return
	
'''
  def load_trajectory(self, trajectory_out_fname):
  
    # We will be loading in the trajectory header first, then 1 frame at a time while self.display_flags['load_trajectory'] == 1
    # Because we are now using pymol, we don't need to store stuff on the blob. Delete from traj as soon as successfully drawn
    
    # Header first
    traj = FFEA_trajectory.FFEA_trajectory(fname = trajectory_out_fname, load_all=0)
  	
  	# Now frames
  	while(True):
  	
  		if self.display_flags['load_trajectory'] == 1:
	  		success = traj.load_frame()
	  		
	  		if success == 0:
	  		
	  			# Success! Copy the frames into the blobs
	  			for i in range(self.script.params.num_blobs):
	  				for j in range(self.script.params.num_conformations[i]):
	  					self.blob_list[i][j].frames.append(traj.blob[i][j].frame[-1])
	  					self.blob_list[i][j].num_frames += 1
	  					
	  			# self.draw_stuff() # we are now drawing every frame after reading it. 
      			if self.num_frames > 1:
        			cmd.mset("1-"+str(self.num_frames))
        		if self.num_frames > 2:
          			cmd.mplay()
	  		else:
	  		
	  			lowest_num_frames = float("inf")
	  			
	  			# Make sure all blobs have same number of frames
	  			for i in range(self.script.params.num_blobs):
	  				for j in range(self.script.params.num_conformations[i]):
	  					if len(traj.blob[i][j].frame) < lowest_num_frames:
	  						lowest_num_frames = len(traj.blob[i][j].frame)
	  			
	  			for i in range(self.script.params.num_blobs):
	  				for j in range(self.script.params.num_conformations[i]):
	  					traj.blob[i][j].frame = traj.blob[i][j].frame[:lowest_num_frames]
	  			
	  			traj.num_frames = lowest_num_frames
	  			
	  			# Wait a bit
'''	  				
'''
  def load_trajectory(self, trajectory_out_fname, ):

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

          # skip asterisk
          line = traj.readline().rstrip()
          if line != "*":
              print "Missing '*' at start of frame", self.num_frames
              print "Instead found '" + line + "'"
              break

          for i in range(self.num_blobs):

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
              if self.blob_list[i][j].get_state() == "STATIC":
                  traj.readline() # skip "Blob x, Conformation y, step z" line
                  traj.readline() # skip "STATIC" line
                  self.blob_list[i][j].load_frame(None)
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
      	
                      for i in range(self.num_blobs):
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

              # plot the last frame 
              for i in range(self.num_blobs):
                # if self.blob_list[i][j].get_state() == "STATIC":
                #   continue
                j = active_conf[i]
                self.blob_list[i][j].draw_frame(-1, self.display_flags)

              if self.springs != None: # and self.display_flags
                self.draw_springs()

              # delete the frames to save memory, but keep frame 0;
              for i in range(self.num_blobs):
                j = active_conf[i]
                if len(self.blob_list[i][j].frames) > 2: 
                  self.blob_list[i][j].frames.pop() # we will only keep the frames
                                             #  if somebody proves that they are useful

          else:
              break

      traj.close()

      # print "read ", self.num_frames, " frames"
      print "finished loading frames"

      # self.draw_stuff() # we are now drawing every frame after reading it. 
      if self.num_frames > 1:
        cmd.mset("1-"+str(self.num_frames))
        if self.num_frames > 2:
          cmd.mplay()
'''

  def get_system_dimensions(self):

      dims = [[float("inf"), -1 * float("inf")] for i in range(3)]

      for b in self.blob_list:
          bdims = b[0].get_dimensions()

          for i in range(3):
              if bdims[i][0] < dims[i][0]:
                  dims[i][0] = bdims[i][0]
              if bdims[i][1] > dims[i][1]:
                  dims[i][1] = bdims[i][1]
      	
      return dims



  def init_vars(self):
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

        self.display_flags = {'show_mesh': 0, ## PYMOL OK
            'show_solid': 1, ## PYMOL OK
            'show_flat': 0,
            'show_material': 0,
            'show_vdw_only': 0,
            'show_node_numbers': 0, ## PYMOL OK 
            'show_pinned_nodes': 0,
            'hide_frozen': 0,
            'show_shortest_edge': 0,
            'vdw_edit_mode': 0,
            'binding_site_edit_mode': 0,
            'binding_site_list': [],
            'selected_index': 0,
            'selected_blob': 0,
            'selected_conformation':0,
            'show_linear_nodes_only': 0,
            'show_mesh_surf': 0, ## PYMOL OK
            'load_trajectory': 1, ## PYMOL OK
            'show_inverted': 0,
            'blob_colour': (1.0, 1.0, 1.0)}

        self.buttons = {'show_mesh' : self.show_mesh,
            'show_solid' : self.show_solid,
            'show_mesh_surf' : self.show_mesh_surf}

        # start the buttons to the default values
        for k in self.buttons.keys():
          self.buttons[k].set(self.display_flags[k])

        self.selected_index = 0
        self.selected_blob = 0
        self.selected_conformation = 0

        self.offset_x = 0
        self.offset_y = 0
        self.offset_z = 0

        self.box_x = -1
        self.box_y = -1
        self.box_z = -1
        self.springs = None
        self.show_box = 0;

        self.modifying_frame = False

        self.recording = 0
        self.movie_dir = "__temp__FFEA_viewer_movie_dir__"

        self.projection = "perspective"



  def draw_stuff(self):

    for f in range(self.num_frames):
      for i in range(self.num_blobs):
          for j in range(self.num_conformations[i]):
              self.blob_list[i][j].draw_frame(f, self.display_flags) ## PLUGIN OUT


      if self.springs != None: # and self.display_flags
         self.draw_springs()


  def draw_springs(self):

      for s in self.springs.spring:

         # Get correct frames
         correct_frame = [-1 for i in range(self.num_blobs)]
         for i in range(self.num_blobs):
            if self.blob_list[i][0].state == "STATIC":
               correct_frame[i] = 0
         print "correct_frame: ", correct_frame

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
         


