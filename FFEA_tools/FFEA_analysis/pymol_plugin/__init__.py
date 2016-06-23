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
     self.root.geometry("155x180")
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
     check_button_show_springs.select() # that has to match with the default value 1! 
     
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
     # WARNINGs:
     NOT_IMPLEMENTED = ["show_element_numbers", "show_face_numbers"]
     if NOT_IMPLEMENTED.count(key):
       print key, " functionality is still under development."
       


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
  # we will take the comments out of iFile,
  #   write a virtual file "ffea_in" 
  #   and return its handler.
  # # # # # # # # # # # # # # # # # # # # # #
  def commentsOut(self, iFile):
     sta = open(iFile, 'r')
     STA = sta.readlines()
     sta.close()

     ffea_in = StringIO.StringIO()

     # and some variables to take the comments out: 
     comment = 0
     m_ini = "<!--"
     m_end = "-->"
     # Now start parsing the input file
     for txt in STA:
     
         # Strip tag wrapping 
         # line = ffea_in.readline().strip()
         line = txt.strip()
     
         # The following stuff takes care of the comments enclosed in "<!--" and "-->":
         # buf2_string = ""
         found = 0
         count = 0
         count_0 = 0
         ini = 0
         end = len(line)
         theEnd = end
     
         # remove the comments:
         while ((found != -1) and (found != len(line))):
           if (comment == 0):
             found = line.find(m_ini)
             if (found != -1):
               count += 1
               comment = 1
               ini = found
           if (comment == 1):
             found = line.find(m_end)
             if (found != -1):
               count += 2
               comment = 0
               end = found + 3
           # the line end up without closing the comment:
           if (comment == 1):
             line = line[:ini]
             break
           # we're out of the comment:
           elif (comment == 0):
             if (count == count_0 + 3):
               buf2_string = line[:ini]
               buf2_string += line[end:]
               line = buf2_string
               count_0 = count
             elif (count == count_0 + 2):
               line = line[end:]
               count_0 = count
         # comments removed!
         if len(line) > 0:
           ffea_in.write(line.strip() + "\n")

     ffea_in.seek(0,0)   
     return ffea_in



  # # # # # # # # # # # # # # # # # # # # # #
  # # # # # # Load the FFEA file # # # # # # 
  # # # # # # # # # # # # # # # # # # # # # # 
  def load_ffea(self, ffea_fname):
     print "diplay flags:", self.display_flags
     print "show springs: ", self.display_flags['show_springs']    

     # Check if given file exists
     if os.path.isfile(ffea_fname) == False:
             print "No such file:", ffea_fname
             return
     else: 
        self.ffea_fname = ffea_fname

     print "Loading ffea file: " + self.ffea_fname
     ffea_path, ffea_id_string = os.path.split(self.ffea_fname)
     if ffea_path == "":
         ffea_path = "."

     # First, we want to see if there is a parameters file associated with this
     # ffea_in = open(self.ffea_fname, "r")
     ffea_in = self.commentsOut(self.ffea_fname)
     	
     # Read required stuff from params block
     trajectory_out_fname = None
     kappa = None
     es_N_x = None
     es_N_y = None
     es_N_z = None
     es_h = None
     self.calc_vdw = 0
     self.move_into_box = 1
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

             ### PLUGIN OUT # not loading the trajectory!
             if self.display_flags['load_trajectory'] == 0:
               print "requested not to load the trajectory"
               trajectory_out_fname = None

         elif lvalue == "calc_vdw":
             self.calc_vdw = int(rvalue)
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
             print "read num_blobs:", self.num_blobs
         elif lvalue == "move_into_box":
             self.move_into_box = int(rvalue)
     			
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
             print "num_conformations: ", self.num_conformations
             print "num_blobs: ", self.num_blobs
             if len(self.num_conformations) != self.num_blobs:
                 sys.exit("Error. Not enough specified 'num_conformations' to create blob array.")

             # Get blob list array
             self.blob_list = [[None for j in range(self.num_conformations[i])] for i in range(self.num_blobs)]

     self.box_x = (1.0 / kappa) * es_h * es_N_x
     self.box_y = (1.0 / kappa) * es_h * es_N_y
     self.box_z = (1.0 / kappa) * es_h * es_N_z
     

     # Let control window know about num_blobs and num_conformations ## PLUGIN OUT
     # self.speak_to_control.send({'num_blobs': self.num_blobs, 'num_conformations': self.num_conformations, 'num_frames': self.num_frames, 'current_frame': -1, 'death': False, 'pausing':self.pausing})

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
#          elif line == "spring":
#              # Ignore
#              while True:
#                  line = ffea_in.readline().strip()[1:-1]
#                  if line == "/spring":
#                      break
#                  else:
#                      continue
         elif line == "interactions":
            while True:
               line = ffea_in.readline().strip()[1:-1]

               if line == "/interactions":
                  print line
                  break
               elif line == "springs" or line == "spring":

                  while True:
                     line = ffea_in.readline().strip()[1:-1]
                     print line
                     if line == "/springs" or line == "/spring":
                        break
                     else:
                        # Ignore
                        #continue

                        # Don't ignore
                        lvalue, rvalue = line.split("=")
                        lvalue = lvalue.strip()
                        rvalue = rvalue.strip()
                        if lvalue == "springs_fname" or lvalue == "spring_fname":
                           if os.path.isabs(rvalue) == False:
                                                rvalue = os.path.join(ffea_path, rvalue)

                           self.springs = FFEA_springs.FFEA_springs(rvalue)
                           self.springs.print_details()

               else: # Add beads block here maybe?
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
                 # new_blob = Blob.Blob(energy_thresh=self.energy_threshold)
                 new_blob = Blob.Blob()
                 new_blob.load(blob_number, blob_index, conformation_index, blob_nodes[i], blob_top[i], blob_surface[i], blob_vdw[i], scale, blob_motion_state[i], blob_pin[i], blob_mat[i], blob_binding[i], blob_centroid_pos, blob_rotation, ffea_path)

                 self.blob_list[blob_index][i] = new_blob
                 new_blob_name = ffea_id_string + "#" + str(blob_index) + ", " + str(conformation_index)
                 info_string = "Name:\t" + ffea_id_string + "\nConformation:\t" + str(i) + "\nNodes:\t" + blob_nodes[i] + "\nTopology:\t" + str(blob_top[i]) + "\nSurface:\t" + blob_surface[i] + "\nVdW:\t" + str(blob_vdw[i]) + "\npin:\t" + str(blob_pin[i]) + "\nMotion State:\t" + blob_motion_state[i] + "\n"
                 add_blob_info = {'name': new_blob_name, 'info': info_string}
                 # self.speak_to_control.send({'add_blob': add_blob_info}) ## PLUGIN OUT 
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

     # Make everything have unitish sizes
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

         b[0].load_nodes_file_as_frame()
	 print b[0].frames[0].node_list[0]
         x, y, z = b[0].get_centroid(0)
         world_centroid[0] += x * b[0].num_nodes
         world_centroid[1] += y * b[0].num_nodes
         world_centroid[2] += z * b[0].num_nodes
         total_num_nodes += b[0].num_nodes

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
         elif b[0].state == "STATIC":
             if self.calc_vdw == 1 and self.move_into_box == 1:
                 b[0].frames[0].translate(shift)
                 print "--- translated by: ", shift
         else:
             b[0].frames = []
             b[0].num_frames = 0

     ## PLUGIN OUT # the following two lines are out for the moment. 
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

              if self.springs != None and (self.display_flags['show_springs']): # and self.display_flags
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
            'show_springs':1, ## PYMOL work in progress
            'show_face_numbers':0, ## PYMOL stub
            'show_element_numbers':0, # PYMOL stub
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


      if (self.springs != None) and (self.display_flags['show_springs']): # and self.display_flags
         self.draw_springs()


  def draw_springs(self):

      for s in self.springs.spring:

         # Get correct frames
         correct_frame = [-1 for i in range(self.num_blobs)]
         for i in range(self.num_blobs):
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
         


