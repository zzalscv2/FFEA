# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#
import cProfile

import sys, os, time
import numpy as np

from pymol import cmd
if (cmd.get_version()[1] < 1.7):
   print "You are running PyMOL v ", cmd.get_version()[1], " but the FFEAPlugin needs a version higher than 1.7"
   raise Exception("FFEAPlugin --- You need a version higher than 1.7 to run the FFEAPlugin")

from pymol.callback import Callback
import warnings

try:
    from mtTkinter import *
except ImportError:
    warnings.warn("DANGER: Tkinter is not thread-safe. PyMOL will eventually crash while loading an FFEA trajectory. Please install mtTKinter.", RuntimeWarning)
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
from pymol import CmdException, stored
from pymol.cgo import *
from pymol.vfont import plain

# FFEA stuff
import FFEA_script
import FFEA_trajectory
import FFEA_turbotrajectory
import FFEA_surface
import FFEA_pin, FFEA_vdw, FFEA_lj
import FFEA_rod

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
     
     Pmw.initialise(self.root)

     self.root.geometry("500x400")
     
     self.root.title("FFEA Loader")

     top_frame = Frame(self.root)
     top_frame.pack()
     


#     menubar = Menu(top_frame)

#     filemenu = Menu(menubar, tearoff=0)
#     filemenu.add_command(label="Load 'ffea' file", command=self.choose_ffea_file_to_load)
#     menubar.add_cascade(label="File", menu=filemenu)
#     self.root.config(menu=menubar)

     # PLUGIN (separated into mutually exclusive sets. Devs take note!)



     ## devs, please use init_vars to centralise the initialisation of values
     self.init_vars()
     self.system_name = StringVar(self.root, value=self.display_flags['system_name'])
     self.do_load_trajectory = StringVar(self.root, value=self.display_flags['load_trajectory'])
     self.show_box = StringVar(self.root, value=self.display_flags['show_box'])
     self.show_pinned = IntVar(self.root, value=self.display_flags['show_pinned'])
     self.show_beads = StringVar(self.root, value=self.display_flags['show_beads'])
     self.show_danger = IntVar(self.root, value=self.display_flags['show_danger'])
     self.show_inverted = IntVar(self.root, value=self.display_flags['show_inverted'])
     self.show_springs = IntVar(self.root, value=self.display_flags['show_springs'])
     self.show_numbers = StringVar(self.root, value=self.display_flags['show_numbers'])
     self.matparam = StringVar(self.root, value=self.display_flags['matparam'])
     self.show_mesh = StringVar(self.root, value=self.display_flags['show_mesh'])
     self.show_shortest_edge = IntVar(self.root, value=self.display_flags['show_shortest_edge'])
     self.load_sfa = StringVar(self.root, value=self.display_flags['load_sfa'])
     self.highlight = StringVar(self.root, value=self.display_flags['highlight'])

     self.sele_name = StringVar(self.root, value=self.display_flags['sele_name'])
     self.pin_fname = StringVar(self.root, value=self.display_flags['pin_fname'])
     self.mat_fname = StringVar(self.root, value=self.display_flags['mat_fname'])
     self.which_vdw_type = IntVar(self.root, value=self.display_flags['which_vdw_type'])
     self.vdw_type0 = IntVar(self.root, value=self.display_flags['vdw_type0'])
     self.vdw_type1 = IntVar(self.root, value=self.display_flags['vdw_type1'])
     self.vdw_fname = StringVar(self.root, value=self.display_flags['vdw_fname'])
     self.lj_fname = StringVar(self.root, value=self.display_flags['lj_fname'])
     self.lj_eps = StringVar(self.root, value=str(self.display_flags['lj_eps']))
     self.lj_r0 = StringVar(self.root, value=str(self.display_flags['lj_r0']))
     self.mat_d = StringVar(self.root, value=str(self.display_flags['mat_d']))
     self.mat_sm = StringVar(self.root, value=str(self.display_flags['mat_sm']))
     self.mat_bm = StringVar(self.root, value=str(self.display_flags['mat_bm']))
     self.mat_sv = StringVar(self.root, value=str(self.display_flags['mat_sv']))
     self.mat_bv = StringVar(self.root, value=str(self.display_flags['mat_bv']))
 
     # # Display flags frame
     
     self.notebook = Pmw.NoteBook(self.root)
     self.notebook.pack(fill = 'both', expand = 1, padx = 10, pady = 10)
     page = self.notebook.add('Loader')
     self.notebook.tab('Loader').focus_set()

     
     display_flags_frame = Frame(page)
     display_flags_frame.pack(anchor=CENTER, expand=True)


     # propose a system name:
     label_system_name = Label(display_flags_frame, text="System name:")
     label_system_name.grid(row=0, column=0, sticky=E)
     self.text_button_system_name = Entry(display_flags_frame, text="load as:", textvariable=self.system_name, validate="focus", validatecommand=lambda:self.update_display_flags("system_name", val=-2, text=self.system_name.get()))
     self.text_button_system_name.grid(row=0, column=1, sticky=W)
     
     self.random_name_button = Button(display_flags_frame, text="Random Name", command=lambda:self.new_system_name());
     self.random_name_button.grid(row=0, column=2, sticky=W)

     label_display = Label(display_flags_frame, text="Display:")
     label_display.grid(row=2, column=2, sticky=W)

     # show springs: 
     self.check_button_show_springs = Checkbutton(display_flags_frame, text="Springs", variable=self.show_springs, command=lambda:self.update_display_flags("show_springs"))
     self.check_button_show_springs.grid(row=3, column=2, sticky=W)


     # show pinned_nodes: 
     self.check_button_show_pinned = Checkbutton(display_flags_frame, text="Pinned Nodes", variable=self.show_pinned, command=lambda:self.update_display_flags("show_pinned"))
     self.check_button_show_pinned.grid(row=4, column=2, sticky=W)


     # show inverted_elements: 
     self.check_button_show_inverted = Checkbutton(display_flags_frame, text="Inverted Elements", variable=self.show_inverted, command=lambda:self.update_display_flags("show_inverted"))
     self.check_button_show_inverted.grid(row=5, column=2, sticky=W)

     # show danger_elements: 
     self.check_button_show_danger = Checkbutton(display_flags_frame, text="Dangerous Elements", variable=self.show_danger, command=lambda:self.update_display_flags("show_danger"))
     self.check_button_show_danger.grid(row=6, column=2, sticky=W)

     # # show solid:
     label_solid = Label(display_flags_frame, text="Show Solid:")
     label_solid.grid(row=1, column=0, sticky=E)
     # Selectable box for material param, i. e., show solid:
     self.spinbox_material_param = OptionMenu(display_flags_frame, self.matparam, "Plain Solid", "Density", "Shear Viscosity", "Bulk Viscosity", "Shear Modulus", "Bulk Modulus", "VdW", "No Solid", command=lambda x:self.update_display_flags("matparam", val=self.matparam.get()) )
     self.spinbox_material_param.grid(row=1, column=1, sticky=W)


     # # show mesh:
     label_mesh = Label(display_flags_frame, text="Show Mesh:")
     label_mesh.grid(row=2, column=0, sticky=E)
     self.om_show_mesh = OptionMenu(display_flags_frame, self.show_mesh, "Surface Mesh", "Whole Mesh", "No Mesh", command=lambda x:self.update_display_flags("show_mesh", val=self.show_mesh.get()))
     self.om_show_mesh.grid(row=2, column=1, sticky=W)


     # # show Numbers:
     label_mesh = Label(display_flags_frame, text="Show Indices:")
     label_mesh.grid(row=3, column=0, sticky=E)
     self.index_option = OptionMenu(display_flags_frame, self.show_numbers, "Node Indices", "Node Indices (Linear)", "Element Indices", "Face Indices", "No Indices", command=lambda x:self.update_display_flags("show_numbers", val=self.show_numbers.get()) )
     self.index_option.grid(row=3, column=1, sticky=W)
     
    
     # Outer simulation box
     label_box = Label(display_flags_frame, text="Show Box:")
     label_box.grid(row=4, column=0, sticky=E)
     self.om_show_box = OptionMenu(display_flags_frame, self.show_box, "Simulation Box (outline)", "Simulation Box (whole)", "No Box", command=lambda x:self.update_display_flags("show_box", val=self.show_box.get()))
     self.om_show_box.grid(row=4, column=1, sticky=W)


     # CG Beads:
     label_beads = Label(display_flags_frame, text="Show Beads:")
     label_beads.grid(row=5, column=0, sticky=E)
     self.om_show_beads = OptionMenu(display_flags_frame, self.show_beads, "Configuration", "Configuration & Assignments", "Trajectory", "No Beads", command=lambda x:self.update_display_flags("show_beads", val=self.show_beads.get()))
     self.om_show_beads.grid(row=5, column=1, sticky=W)


     ## # Trajectory Dropdown Menu # #
     label_traj = Label(display_flags_frame, text="Load:")
     label_traj.grid(row=6, column=0, sticky=E)
     self.om_do_load_trajectory = OptionMenu(display_flags_frame, self.do_load_trajectory, "Trajectory", "System (Into box)", "System (Plainly)", "CGO", command=lambda x:self.update_display_flags("load_trajectory", val=self.do_load_trajectory.get())) 
     self.om_do_load_trajectory.grid(row=6, column=1, sticky=W)


     ## # Add Supportive Fake Atoms (SFA) box # #
     label_sfa = Label(display_flags_frame, text="Add Atoms:")
     label_sfa.grid(row=7, column=0, sticky=E)
     self.om_load_sfa = OptionMenu(display_flags_frame, self.load_sfa, "None", "Onto Linear Nodes", "Onto Nodes", "Onto Faces", "Onto Elements", command=lambda x:self.update_display_flags("load_sfa", val=self.load_sfa.get())) 
     self.om_load_sfa.grid(row=7, column=1, sticky=W)
     
     
     ## # Finally the Load Button! # #
     self.load_button = Button(display_flags_frame, text="Load ffea file", command=lambda:self.choose_ffea_file_to_load() )
     self.load_button.grid(row=8, column=0, columnspan=4, sticky=W+E+N+S, pady=20)
     



     # # # # EDITOR TAB # # # # # 
     page = self.notebook.add('Editor')
     self.notebook.tab('Editor').focus_set()
     
     editor_frame = Frame(page)
     editor_frame.pack(anchor=CENTER, fill='both',expand=True)

     ## ## ## Selection Frame and Box ## ## ##
     sele_frame = Frame(editor_frame)
     sele_frame.pack(fill=X, side=TOP)

     sele_label = Label(sele_frame, text="Selection:")
     sele_label.pack(side=LEFT)
     self.text_button_sele_name = Entry(sele_frame, text="sele", textvariable=self.sele_name, validate="focus", validatecommand=lambda:self.update_display_flags("sele_name", val=-2, text=self.sele_name.get()))
     self.text_button_sele_name.pack(side=LEFT)
     self.text_button_sele_name.config(state=DISABLED)
    
     ## ## ## Material Box ## ## ##
     mat_group = Pmw.Group(editor_frame,tag_text='Material Parameters')
     mat_group.pack(fill='both' ,side=TOP)

     ## Contents
     # Mat Params
     label_mat_d = Label(mat_group.interior(), text="Density:")
     label_mat_d.grid(row=0,column=0)
     self.text_button_mat_d = Entry(mat_group.interior(), width=12, text="1.5e3", textvariable=self.mat_d, validate="focus", validatecommand=lambda:self.update_display_flags("mat_d", val=-3, text=self.mat_d.get()))
     self.text_button_mat_d.grid(row=0,column=1)
     self.text_button_mat_d.config(state=DISABLED)
     label_mat_sm = Label(mat_group.interior(), text="Shear Modulus:")
     label_mat_sm.grid(row=1,column=0)
     self.text_button_mat_sm = Entry(mat_group.interior(), width=12, text="370.37e6", textvariable=self.mat_sm, validate="focus", validatecommand=lambda:self.update_display_flags("mat_sm", val=-3, text=self.mat_sm.get()))
     self.text_button_mat_sm.grid(row=1,column=1)
     self.text_button_mat_sm.config(state=DISABLED)
     label_mat_bm = Label(mat_group.interior(), text="Bulk Modulus:")
     label_mat_bm.grid(row=2,column=0)
     self.text_button_mat_bm = Entry(mat_group.interior(), width=12, text="111.11e7", textvariable=self.mat_bm, validate="focus", validatecommand=lambda:self.update_display_flags("mat_bm", val=-3, text=self.mat_bm.get()))
     self.text_button_mat_bm.grid(row=2,column=1)
     self.text_button_mat_bm.config(state=DISABLED)
     label_mat_sv = Label(mat_group.interior(), text="Shear Viscosity:")
     label_mat_sv.grid(row=3,column=0)
     self.text_button_mat_sv = Entry(mat_group.interior(), width=12, text="1e-3", textvariable=self.mat_sv, validate="focus", validatecommand=lambda:self.update_display_flags("mat_sv", val=-3, text=self.mat_sv.get()))
     self.text_button_mat_sv.grid(row=3,column=1)
     self.text_button_mat_sv.config(state=DISABLED)
     label_mat_bv = Label(mat_group.interior(), text="Bulk Viscosity:")
     label_mat_bv.grid(row=4,column=0)
     self.text_button_mat_bv = Entry(mat_group.interior(), width=12, text="1e-3", textvariable=self.mat_bv, validate="focus", validatecommand=lambda:self.update_display_flags("mat_bv", val=-3, text=self.mat_bv.get()))
     self.text_button_mat_bv.grid(row=4,column=1)
     self.text_button_mat_bv.config(state=DISABLED)

     self.load_mat_button = Button(mat_group.interior(), text="Update Material file", command=lambda:self.choose_mat_file_to_setup() )
     self.load_mat_button.grid(column=3, row=2)
     self.load_mat_button.config(state=DISABLED)   

     ## ## ## Pin Box ## ## ## 
     pin_group = Pmw.Group(editor_frame,tag_text='Pinned Nodes')
     pin_group.pack(fill='both' ,side=TOP)

     ## Contents
     self.load_pin_button = Button(pin_group.interior(), text="Update Pin file", command=lambda:self.choose_pin_file_to_setup());
     self.load_pin_button.pack(side=LEFT)
     self.load_pin_button.config(state=DISABLED)

     ## ## ## VdW & LJ Box ## ## ##
     vdwlj_group = Pmw.Group(editor_frame,tag_text='Van der Waals Interactions')
     vdwlj_group.pack(fill='both', side=TOP)

     ## Contents
     # VdW Type frame
     vdwtype_frame = Frame(vdwlj_group.interior())
     vdwtype_frame.pack(fill=X, side=TOP)

     # Define vdw face types and a radio button to switch between them:
     self.radio_vdw_type0 = Radiobutton(vdwtype_frame, variable=self.which_vdw_type, value=0, command=lambda:self.update_display_flags("which_vdw_type", val=0))
     self.radio_vdw_type0.pack(side=LEFT, fill=X)
     label_vdw_type0 = Label(vdwtype_frame, text="VdW type 1:")
     label_vdw_type0.pack(side=LEFT)
     self.spinbox_vdw_type0 = OptionMenu(vdwtype_frame, self.vdw_type0, "-1 (no vdw)", "0", "1", "2", "3", "4", "5", "6", "7", command=lambda x: self.update_display_flags("vdw_type0", val=self.vdw_type0.get()))
     self.spinbox_vdw_type0.pack(side=LEFT)
     self.spinbox_vdw_type0.config(state=DISABLED)

     self.spinbox_vdw_type1 = OptionMenu(vdwtype_frame, self.vdw_type1, "-1 (no vdw)", "0", "1", "2", "3", "4", "5", "6", "7", command=lambda x: self.update_display_flags("vdw_type1", val=self.vdw_type1.get()))
     self.spinbox_vdw_type1.pack(side=RIGHT)
     self.spinbox_vdw_type1.config(state=DISABLED)
     label_vdw_type1 = Label(vdwtype_frame, text="VdW type 2:")
     label_vdw_type1.pack(side=RIGHT)
     self.radio_vdw_type1 = Radiobutton(vdwtype_frame, variable=self.which_vdw_type, value=1, command=lambda:self.update_display_flags("which_vdw_type", val=1))
     self.radio_vdw_type1.pack(side=RIGHT, fill=X)

     ## ## ## VdW Box ## ## ##
     vdw_group = Pmw.Group(vdwlj_group.interior(),tag_text='VdW Files')
     vdw_group.pack(fill='both', expand=1, side=LEFT)

     ## Contents
     
     # choose a vdw file to set up 
     self.load_vdw_button = Button(vdw_group.interior(), text="Update VdW file", command=lambda:self.choose_vdw_file_to_setup() )
     self.load_vdw_button.pack(anchor=CENTER)
     self.load_vdw_button.config(state=DISABLED)

     ## ## ## LJ Box ## ## ##
     lj_group = Pmw.Group(vdwlj_group.interior(),tag_text='LJ Files')
     lj_group.pack(fill='both', expand=1, side=LEFT)

     ## Contents
     # LJ Params frame
     ljparams_frame = Frame(lj_group.interior())
     ljparams_frame.pack(fill=X, side=TOP)

     # Define lj params:
     label_lj_eps = Label(ljparams_frame, text="LJ Epsilon:")
     label_lj_eps.pack(side=LEFT, fill=X)
     self.text_button_lj_eps = Entry(ljparams_frame, width=12, text="1e15", textvariable=self.lj_eps, validate="focus", validatecommand=lambda:self.update_display_flags("lj_eps", val=-3, text=self.lj_eps.get()))
     self.text_button_lj_eps.pack(side=LEFT)
     self.text_button_lj_eps.config(state=DISABLED)

     self.text_button_lj_r0 = Entry(ljparams_frame, width=12, text="1e-9", textvariable=self.lj_r0, validate="focus", validatecommand=lambda:self.update_display_flags("lj_r0", val=-3, text=self.lj_r0.get()))
     self.text_button_lj_r0.pack(side=RIGHT)
     self.text_button_lj_r0.config(state=DISABLED)
     label_lj_r0 = Label(ljparams_frame, text="LJ r0:")
     label_lj_r0.pack(side=RIGHT, fill=X)

     self.load_lj_button = Button(lj_group.interior(), text="Update LJ file", command=lambda:self.choose_lj_file_to_setup() )
     self.load_lj_button.pack(side=TOP)
     self.load_lj_button.config(state=DISABLED)

     ## ## ## To be called after using Pmw.Group:
     self.notebook.setnaturalsize()

 
  def update_material(self):

	print("Updating Material File...\n")

	## Check we're ok and ready to go
	# Is there a selection?
	try:
		num_atoms = cmd.count_atoms(self.display_flags["sele_name"])
	except CmdException:
		print("Cannot make material file as '" + self.display_flags["sele_name"] + "' selection does not exist.")
		return

	# Does it have atoms?
	if num_atoms == 0:
		print("Will not make material file as '" + self.display_flags["sele_name"] + "' contains 0 selected pseudoatoms.")
		return

	# Are they all in the same blob?
	stored.blob_IDs = []
	cmd.iterate(self.display_flags["sele_name"], "stored.blob_IDs.append(model)")
	if stored.blob_IDs.count(stored.blob_IDs[0]) != len(stored.blob_IDs):
		print("Cannot make material file as '" + self.display_flags["sele_name"] + "' contains pseudoatoms from more than 1 blob.")
		return

	# Get the blob_ID, to identify the blob within self.blob_list
	blob_ID = int(stored.blob_IDs[0].split("_")[1])

	# If the blob is not DYNAMIC, abort:
	if self.blob_list[blob_ID][0].motion_state != 'DYNAMIC':
		print("Selection :" + self.display_flags["sele_name"] + " belongs to a non-dynamic blob. Because of being static, material parameters would not make any effect, so we are aborting. Edit your ffea file to make the corresponding blob DYNAMIC if you want to pursue.")
		return


	## Ok, we're ready!
	
	# Get indices and index type
	stored.baseindices = []
	cmd.iterate(self.display_flags["sele_name"], "stored.baseindices.append(resi)")
	stored.baseindices = [int(i) for i in stored.baseindices]
	indextype = stored.blob_IDs[0].split("_")[2]

	# Different things depending on index type
	indices = []
	if indextype == "efa":

		# Element indices. Great!
		indices = stored.baseindices

	elif indextype == "nfa" or indextype == "lnfa":

		# Node indices. Get elements with 1 or more indices present
		try:
			indices = self.blob_list[blob_ID][0].top.index_switch(stored.baseindices, "node", limit=1)
		except(IndexError):
			print("Could not make material file as correct indices could not be extracted from topology using node selection")
			return

	elif indextype == "sfa" or indextype == "ffa":

		# Surface indices. Get elements with surface face present
		try:
			indices = self.blob_list[blob_ID][0].top.index_switch(stored.baseindices, "surf", limit=1, surf=self.blob_list[blob_ID][0].surf)
		except(IndexError):
			print("Could not make material file as correct indices could not be extracted from topology using face selection")
			return

		except(IOError):
			print("Could not make material file as there was a problem linking topology to surface")

	# Get a new material file and copy relevent stuff
	mat = FFEA_material.FFEA_material()
	mat.element = self.blob_list[blob_ID][0].mat.element
	mat.num_elements = self.blob_list[blob_ID][0].mat.num_elements

	# Update stuff
	for i in indices:
		mat.set_params(i, self.display_flags["mat_d"], self.display_flags["mat_sv"], self.display_flags["mat_bv"], self.display_flags["mat_sm"], self.display_flags["mat_bm"], 1.0)

	# Write out
	mat.write_to_file(self.display_flags["mat_fname"])
	print("...done!\n")

  def update_pin(self):

	print("Updating Pin File...\n")

	## Check we're ok and ready to go
	# Is there a selection?
	try:
		num_atoms = cmd.count_atoms(self.display_flags["sele_name"])
	except CmdException:
		print("Cannot make pin file as '" + self.display_flags["sele_name"] + "' selection does not exist.")
		return

	# Does it have atoms?
	if num_atoms == 0:
		print("Will not make pin file as '" + self.display_flags["sele_name"] + "' contains 0 selected pseudoatoms.")
		return

	# Are they all in the same blob?
	# # get the list of names:
	stored.blob_IDs = []
	cmd.iterate(self.display_flags["sele_name"], "stored.blob_IDs.append(model)")
	# # ignore those coming from "pinned" selections:
	ndx = len(stored.blob_IDs) - 1
	for i in reversed(stored.blob_IDs):
		if i.split("_")[2] == "pinned": stored.blob_IDs.pop(ndx)
		ndx -= 1
	# # check if they are all equal:
	if stored.blob_IDs.count(stored.blob_IDs[0]) != len(stored.blob_IDs):
		print("Cannot make pin file as '" + self.display_flags["sele_name"] + "' contains pseudoatoms from more than 1 blob.")
		return
	blob_ID = int(stored.blob_IDs[0].split("_")[1])

	# If the blob is not DYNAMIC, abort:
	if self.blob_list[blob_ID][0].motion_state != 'DYNAMIC':
		print("Selection :" + self.display_flags["sele_name"] + " belongs to a non-dynamic blob. Pinning nodes in this case would not make any effect, so we are aborting. Edit your ffea file to make the corresponding blob DYNAMIC if you want to pursue.")
		return


	## Ok, we're ready!
	
	# Get indices and index type
	stored.baseindices = []
	cmd.iterate(self.display_flags["sele_name"], "stored.baseindices.append(resi)")
	stored.baseindices = [int(i) for i in stored.baseindices]
	indextype = stored.blob_IDs[0].split("_")[2]

	# Different things depending on index type
	indices = []
	if indextype == "lnfa":
		# Node indices. Great!
		indices = stored.baseindices

	elif indextype == "nfa":
		# try to remove the non-linear nodes:
		snfa = []
		for i in stored.baseindices:
			# store only linear nodes
			if self.blob_list[blob_ID][0].linear_node_list.count(i):
				indices.append(i)
			else:
				snfa.append(i)
		if len(snfa) > 0:
			print "Only linear nodes can be pinned, but nodes: ", snfa, " are auxilliary, and adding them in the .pin node file has no effect"

	elif indextype == "efa":

		# Element indices. Get all (linear) nodes on each element
		try:
			indices = self.blob_list[blob_ID][0].node.index_switch(stored.baseindices, "topology", top = self.blob_list[blob_ID][0].top)
		except(IndexError):
			print("Could not make pin file as correct indices could not be extracted from nodes using topology selection")
			return

	elif indextype == "sfa" or indextype == "ffa":

		# Surface indices. Get all nodes on each face
		try:
			indices = self.blob_list[blob_ID][0].node.index_switch(stored.baseindices, "surface", surf = self.blob_list[blob_ID][0].surf)
		except(IndexError):
			print("Could not make pin file as correct indices could not be extracted from nodes using surface selection")
			return

	# Add already pinned stuff
	indices.extend(self.blob_list[blob_ID][0].pin.index)

	# Make unique entries
	indices = list(set(indices))

	# Get a new pin file and copy relevent stuff
	pin = FFEA_pin.FFEA_pin()
	for i in indices:
		pin.add_pinned_node(i) 

	# Write out
	pin.write_to_file(self.display_flags["pin_fname"])
	print("...done!\n")

  def update_vdw(self):

	print("Updating VdW File...\n")

	## Check we're ok and ready to go
	# Is there a selection?
	try:
		num_atoms = cmd.count_atoms(self.display_flags["sele_name"])
	except CmdException:
		print("Cannot make VdW file as '" + self.display_flags["sele_name"] + "' selection does not exist.")
		return

	# Does it have atoms?
	if num_atoms == 0:
		print("Will not make VdW file as '" + self.display_flags["sele_name"] + "' contains 0 selected pseudoatoms.")
		return

	# Are they all in the same blob?
	stored.blob_IDs = []
	cmd.iterate(self.display_flags["sele_name"], "stored.blob_IDs.append(model)")
	if stored.blob_IDs.count(stored.blob_IDs[0]) != len(stored.blob_IDs):
		print("Cannot make VdW file as '" + self.display_flags["sele_name"] + "' contains pseudoatoms from more than 1 blob.")
		return
	blob_ID = int(stored.blob_IDs[0].split("_")[1])

	## Ok, we're ready!
	
	# Get indices and index type
	stored.baseindices = []
	cmd.iterate(self.display_flags["sele_name"], "stored.baseindices.append(resi)")
	stored.baseindices = [int(i) for i in stored.baseindices]
	indextype = stored.blob_IDs[0].split("_")[2]

	# Different things depending on index type
	indices = []
	if indextype == "sfa" or indextype == "ffa":

		# Surface indices. Great!
		indices = stored.baseindices

	elif indextype == "nfa" or indextype == "lnfa":

		# Node indices. Get all (linear) nodes on each element
		try:
			indices = self.blob_list[blob_ID][0].surf.index_switch(stored.baseindices, "node", limit=1)
		except(IndexError):
			print("Could not make VdW file as correct indices could not be extracted from surface using node selection")
			return

	elif indextype == "efa":

		# Element indices. Get all (linear) nodes on each element
		try:
			indices = self.blob_list[blob_ID][0].node.index_switch(stored.baseindices, "element", surf = self.blob_list[blob_ID][0].surf)
		except(IndexError):
			print("Could not make VdW file as correct indices could not be extracted from surface using element selection")
			return

	# Get a vdw file
	newvdw = FFEA_vdw.FFEA_vdw()
	newvdw.num_faces = self.blob_list[blob_ID][0].vdw.num_faces
	newvdw.index = self.blob_list[blob_ID][0].vdw.index

	# Update list
	vdwindex = self.display_flags["vdw_type" + str(self.display_flags["which_vdw_type"])]
	for i in indices: 
		newvdw.set_index(int(i), vdwindex)

	# Write
	newvdw.write_to_file(self.display_flags["vdw_fname"])
	print("...done!\n")

  def update_lj(self):

	print("Updating LJ File...\n")

	# No need to check for selections

	# Load an LJ file and copy stuff
	lj = FFEA_lj.FFEA_lj()
	lj.interaction = self.lj.interaction

	# Update stuff
	lj.set_interaction_pair(self.display_flags["vdw_type0"], self.display_flags["vdw_type1"], self.display_flags["lj_eps"], self.display_flags["lj_r0"])

	# Write out
	lj.write_to_file(self.display_flags["lj_fname"])
 	print("...done!\n")

 #################################################
  # # # # Update display_flags from buttons # # #
  # # use val = -3 for floats (Entries)
  # #     val = -2 for strings (Entries)
  # #     val = -1 for binary choices (Checkboxes)
  # #     val > 0, integer, for Radiobuttons 
 #################################################
  def update_display_flags(self, key, val=-1, text=""):
     # If unset (i.e. checkbutton)
     if val == -3:
	try:
	    self.display_flags[key] = float(text)
	    return True
	except(ValueError):
	    return False

     elif val == -2:
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
     print "ffea_fname: ", ffea_fname
     if len(ffea_fname) == 0:
             return

     # load the file
     # cProfile.runctx('self.load_ffea(ffea_fname)', globals(), locals())
     self.load_ffea(ffea_fname)


  def choose_mat_file_to_setup(self):
     # set up the options for the open file dialog box
     options = {}
     options['defaultextension'] = '.mat'
     options['filetypes'] = [('mat files', '.mat'), ('all files', '.*')]
     options['initialdir'] = os.getcwd()
     options['title'] = 'Choose material file'

     # Ask user to select a file
     self.mat_fname.set(tkFileDialog.asksaveasfilename(**options))
     if self.mat_fname.get() == None or self.mat_fname.get().strip() == "":
	return
     self.update_display_flags("mat_fname", -2, self.mat_fname.get())
     self.update_material()

  def choose_pin_file_to_setup(self):
     # set up the options for the open file dialog box
     options = {}
     options['defaultextension'] = '.pin'
     options['filetypes'] = [('pin files', '.pin'), ('all files', '.*')]
     options['initialdir'] = os.getcwd()
     options['title'] = 'Choose pin file'

     # Ask user to select a file
     self.pin_fname.set(tkFileDialog.asksaveasfilename(**options))
     if self.pin_fname.get() == None or self.pin_fname.get().strip() == "":
	return
     self.update_display_flags("pin_fname", -2, self.pin_fname.get())
     self.update_pin()

  def choose_vdw_file_to_setup(self):
     # set up the options for the open file dialog box
     options = {}
     options['defaultextension'] = '.vdw'
     options['filetypes'] = [('vdw files', '.vdw'), ('all files', '.*')]
     options['initialdir'] = os.getcwd()
     options['title'] = 'Choose vdw file'

     # Ask user to select a file
     self.vdw_fname.set(tkFileDialog.asksaveasfilename(**options))
     if self.vdw_fname.get() == None or self.vdw_fname.get().strip() == "":
	return
     self.update_display_flags("vdw_fname", -2, self.vdw_fname.get())
     self.update_vdw()

  def choose_lj_file_to_setup(self):
     # set up the options for the open file dialog box
     options = {}
     options['defaultextension'] = '.lj'
     options['filetypes'] = [('lj files', '.lj'), ('all files', '.*')]
     options['initialdir'] = os.getcwd()
     options['title'] = 'Choose lj file'

     # Ask user to select a file
     self.lj_fname.set(tkFileDialog.asksaveasfilename(**options))

     if self.lj_fname.get() == None or self.lj_fname.get().strip() == "":
	return
     self.update_display_flags("lj_fname", -2, self.lj_fname.get())
     self.update_lj()

  # # # # # # # # # # # # # # # # # # # # # #
  # # # # # # Load the FFEA file # # # # # # 
  # # # # # # # # # # # # # # # # # # # # # # 
  def load_ffea(self, ffea_fname):
      
	tbegin = time.time()
	self.notebook.selectpage("Editor")
  	
	# Update display flags patch (the .get() function got the old spinbox value, so here it's definitely updated)
	self.display_flags['matparam'] = self.matparam.get()

	# Try to reset previous system and update
	self.num_frames = 0

	# Check if given file exists
	if os.path.isfile(ffea_fname) == False:
		print "No such file:", ffea_fname
		return
	else: 
		self.ffea_fname = ffea_fname
        
	print "Loading ffea file: " + self.ffea_fname
    
    # Load script (comments are now removed inside this module, by the way :) )
	self.script = FFEA_script.FFEA_script(self.ffea_fname, fix=True)
	if (self.script.params == None):
		print "Something went wrong reading the FFEA input file", self.ffea_fname
		return
	p = self.script.params
	bl = self.script.blob
	
	# See whether or not to remove traj file (make this better later i.e. rolling loading by storing file pointers)
	if self.display_flags['load_trajectory'] == "System (Into box)" or self.display_flags['load_trajectory'] == "System (Plainly)":
		print "Requested not to load the trajectory"
		#p.trajectory_out_fname = None
	if self.display_flags['load_trajectory'] == "System (Plainly)":
		print "Requested to show the coordinates as they are in the .node(s) file(s)"
		# print "... equivalently, setting < move_into_box = 0 >"  and no PBC:
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
			# try:
			new_blob.load(idnum, bindex, cindex, self.script, self.display_flags)
			# except:
				# print("ERROR: Could not load Blob %d, conformation %d. Please try again." % (bindex, cindex))
				# return

			self.blob_list[bindex][cindex] = new_blob
			new_blob_name = ffea_id_string + "#" + str(bindex) + ", " + str(cindex)
			info_string = "Name:\t" + ffea_id_string + "\nConformation:\t" + str(cindex) + "\nNodes:\t" + c.nodes + "\nTopology:\t" + c.topology + "\nSurface:\t" + c.surface + "\nVdW:\t" + c.vdw + "\npin:\t" + c.pin + "\nMotion State:\t" + c.motion_state + "\n"
			add_blob_info = {'name': new_blob_name, 'info': info_string}
			
			idnum += 1
                 

	# Load some lj
	try:
		self.lj = self.script.load_lj()
	except:
		print("\nERROR: '" + self.script.params.vdw_forcefield_params + "' could not be loaded.")
		print("\nERROR: Could not load system. Please try again.")
		return

	if (not self.lj.valid): 
		print('Something went wrong initialising lennard-jones (lj) parameters')
		print("\nERROR: Could not load system. Please try again.")
		return

	# Load some springs
	if self.display_flags['show_springs'] == 1:
		try:
			self.springs = FFEA_springs.FFEA_springs(self.script.spring)
		except:
			sys.stdout.write("Springs could not be loaded. Continuing...")
			sys.stdout.flush()
			self.springs = None
			
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
            
    # Rescale rods
	if len(self.script.rod) > 0:
		for rod_num in range(len(self.script.rod)):
			self.script.rod[rod_num].scale(self.global_scale)

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

	# Build the box:
	# Do we need to calculate the size of the box? Double the rounded up size of the system
	for i in range(3):
		if p.es_N[i] < 1:
			dims = self.get_system_dimensions(0)
			for j in range(3):
				p.es_N[j] = 2 * int(np.ceil(dims[j] / (self.global_scale*p.vdw_cutoff)))
			break

	self.box = p.vdw_cutoff * p.es_N
	self.box_exists = True
	
		
	# Rescale box
	self.box *= self.global_scale

	# Shift all blobs to center of box if necessary
	shift = 0.5 * self.box - world_centroid
	# if p.calc_vdw == 1 and p.move_into_box == 1:
	if p.move_into_box == 1:
		for b in self.blob_list:
			b[0].frames[0].translate(shift)
			if b[0].beads.pdb != None: b[0].beads.pdb.translate(shift) # beads only work for conf 0
       # move rods into box
    	if len(self.script.rod) > 0:
    		for rod_num in range(len(self.script.rod)):
    			self.script.rod[rod_num].translate(shift)
    		

	# Now, apply PBC if necessary
	# if p.calc_vdw == 1 and self.display_flags['load_trajectory'] != "System (Plainly)":
	if self.display_flags['load_trajectory'] != "System (Plainly)":
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
	if self.display_flags['load_trajectory'] == "CGO":
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
		waitForTrajToLoad = True
		# self.load_trajectory(p.trajectory_out_fname) # serial


	#
	# Print info for the user that won't be deleted from the command line by the trajectory loading
	#
	if self.display_flags['show_springs'] == 1 and self.springs != None:
		if self.script.params.calc_springs == 0 and self.springs.get_num_springs() > 0:
			for b in self.script.blob:
				if b.solver == "CG_nomass":
					print "INFO: Springs have been drawn but calc_springs == 0 in your script. Please change for ffea simulation if you want to use springs."
					break
				
	if waitForTrajToLoad: 
		self.load_trajectory_thread.join()

	# Requires knowledge of whole trajectory
	if self.traj != None and self.display_flags['load_trajectory'] == "Trajectory" and self.wontLoadTraj != 1:
		if self.display_flags["show_inverted"] == 1:  # Show inverted elements
			self.draw_inverted_elements()
		if self.display_flags["show_beads"] == "Trajectory": # Load the trajectory of the beads.
			if self.script.params.beads_out_fname != "":
				beads_traj_name = self.display_flags['system_name'] + "_b"
				cmd.load(self.script.params.beads_out_fname, beads_traj_name)
				cmd.hide("everything", beads_traj_name)
				cmd.show("spheres", beads_traj_name)
			else:
				print "Beads trajectory won't load: beads_out_fname was not defined in the .ffea input file"

	# Load up ye rods
	if len(self.script.rod) > 0:
		for rod_num in range(len(self.script.rod)):
			self.load_rod(self.script.rod[rod_num], rod_num)

   	# Center everything, zoom and sort clipping plane
	cmd.center()
 	cmd.zoom()

	# deactivate load options:
	self.check_button_show_springs.config(state=DISABLED)
	self.check_button_show_pinned.config(state=DISABLED)
	self.check_button_show_inverted.config(state=DISABLED)
	self.check_button_show_danger.config(state=DISABLED)

	self.text_button_system_name.config(state=DISABLED)
	self.random_name_button.config(state=DISABLED)
	self.spinbox_material_param.config(state=DISABLED)
	self.om_show_mesh.config(state=DISABLED)
	self.index_option.config(state=DISABLED)
	self.om_show_box.config(state=DISABLED)
	self.om_load_sfa.config(state=DISABLED)
	self.om_do_load_trajectory.config(state=DISABLED)
	self.load_button.config(state=DISABLED)

	if self.load_sfa.get() != "None":

		# activate Editor options:
		self.text_button_sele_name.config(state="normal")

		self.text_button_mat_d.config(state="normal")
		self.text_button_mat_sm.config(state="normal")
		self.text_button_mat_bm.config(state="normal")
		self.text_button_mat_sv.config(state="normal")
		self.text_button_mat_bv.config(state="normal")
		self.load_mat_button.config(state="normal")

		self.load_pin_button.config(state="normal")

		self.spinbox_vdw_type0.config(state="normal")
		self.spinbox_vdw_type1.config(state="normal")

		self.load_vdw_button.config(state="normal")

		self.text_button_lj_eps.config(state="normal")
		self.text_button_lj_r0.config(state="normal")

		self.load_lj_button.config(state="normal")
	else:
		self.root.destroy()

	print "System loaded in ", time.time() - tbegin, "s."


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

  def load_rod(self, rod, rod_num=0):
    
    def get_avg_lengths(rod):
        # Get scale factor for ei
        avg_pi = 0
        pi = rod.get_p_i(rod.current_r)
        for frame in pi:
            for p in range(len(frame)):
                avg_pi += np.linalg.norm(frame[p])
        avg_pi = avg_pi / (( len(pi)*len(pi[0]) ))
        
        # Get scale factor for m
        avg_m = 0
        for frame in rod.current_m:
            for m in range(len(frame)-1):
                avg_m += np.linalg.norm(frame[m])
        avg_m = avg_m / (( len(rod.current_m)*len(rod.current_m[0]) ))
        
        return avg_m, avg_pi
    
    avg_m, avg_p = get_avg_lengths(rod)
    
    # Scale the material frame to be a similar size to the rod elements
    rod.current_m /= (avg_m/avg_p)/np.sqrt(2)
    
    # units note: radii are arbitrary so far. the *10**10 is to go from SI to angstroms (I should remove this after I add proper scaling)
    for i in range(len(rod.current_r)):
      line = []
      for j in range(len(rod.current_r[i])-1):
        line = line + [9.0, rod.current_r[i][j][0], rod.current_r[i][j][1], rod.current_r[i][j][2], rod.current_r[i][j+1][0], rod.current_r[i][j+1][1], rod.current_r[i][j+1][2], 10, 0, 1, 0, 0, 1, 0 ]
        # material frame in center of each element
        mid_x, mid_y, mid_z = (rod.current_r[i][j][0]+rod.current_r[i][j+1][0])/2, (rod.current_r[i][j][1]+rod.current_r[i][j+1][1])/2, (rod.current_r[i][j][2]+rod.current_r[i][j+1][2])/2
        line = line + [9.0, mid_x, mid_y, mid_z, mid_x+rod.current_m[i][j][0], mid_y+rod.current_m[i][j][1], mid_z+rod.current_m[i][j][2], 5, 0, 0, 1, 0, 0, 1 ]
#        print("frame "+str(i)+" element "+str(j)+"= "+str(line))
#        print(line)
      cmd.load_cgo(line, self.display_flags['system_name']+"_rod_"+str(rod_num), i)


  def call_add_node_pseudoatoms(self):
     if self.display_flags['load_trajectory'] == "System (Plainly)" or self.display_flags['load_trajectory'] == "CGO" or self.wontLoadTraj == 1:
        self.add_node_pseudoatoms_from_nodes()
     elif self.display_flags['load_trajectory'] == "Trajectory":
        self.add_node_pseudoatoms()
     else:
        print "Cannot add pseudoatoms if selecting System (Into box)" 

    
  def add_node_pseudoatoms(self):
      node_object_list = []
      for blob_num in range(len(self.blob_list)):
          traj = self.script.load_trajectory(1)
          node_object_list.append(traj.blob[blob_num])
      for node_object in range(len(node_object_list)):
          for conformation in node_object_list[node_object]:
              for node in range(len(conformation.frame[0].pos)):
                  if conformation.frame[0].pos[node] !=None:
                      cmd.pseudoatom(pos = (conformation.frame[0].pos[node]*1000000000).tolist(), name = str(node), color="black")
                          
  def add_node_pseudoatoms_from_nodes(self):
      node_object_list = []
      for blob_num in range(len(self.blob_list)):
          node_object_list.append(self.script.load_node(blob_num))
      for node_object in range(len(node_object_list)):
          for node in range(len(node_object_list[node_object].pos)):
              cmd.pseudoatom(pos = node_object_list[node_object].pos[node].tolist())   
        
 
  def draw_inverted_elements(self):

	# For each blob
	bin = 0
	for b in self.blob_list:
		# Change when conformations are stable
		cin = 0
		c = b[cin]

		element_list = []
		
		# Get last two frames and check whether volume / jacobian has changed it's sign
		
		index = 0
		if (c.top == None):
			if (c.motion_state != "STATIC"):
				print("Cannot draw inverted elements for blob %d as there is not topology" % (bin))
			bin += 1
			continue

		flast = self.traj.blob[bin][cin].frame[-1]

		try:
			f2last = self.traj.blob[bin][cin].frame[-2]
		except:
			f2last = c.node

		for el in c.top.element:
			jac = np.linalg.det(el.calc_jacobian(flast))
			jac_last = np.linalg.det(el.calc_jacobian(f2last))
			if jac * jac_last < 0:
				element_list.append(index)

			index += 1
		
		# Draw these as a new object on the last frame		
		invele = []
		numtxt = []
		txtscale = 0.1
		axes = np.array([[15.0,0.0,0.0],[0.0,15.0,0.0],[0.0,0.0,15.0]])
		invele.extend( [BEGIN, LINES] )

		for el in element_list:
			n1 = flast.pos[c.top.element[el].n[0]]
			n2 = flast.pos[c.top.element[el].n[1]]
			n3 = flast.pos[c.top.element[el].n[2]]
			n4 = flast.pos[c.top.element[el].n[3]]

			invele.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
			invele.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )

			invele.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
			invele.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )

			invele.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )
		        invele.extend( [ VERTEX, n4[0], n4[1], n4[2] ] )

		        invele.extend( [ VERTEX, n4[0], n4[1], n4[2] ] )
		        invele.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )

		        invele.extend( [ VERTEX, n1[0], n1[1], n1[2] ] )
			invele.extend( [ VERTEX, n3[0], n3[1], n3[2] ] )

			invele.extend( [ VERTEX, n2[0], n2[1], n2[2] ] )
		        invele.extend( [ VERTEX, n4[0], n4[1], n4[2] ] )

			nn = c.top.element[el].calc_centroid(flast)
			cyl_text(numtxt,plain,nn,str(el), txtscale, axes=axes * txtscale)

		invele.append(END)

		if len(invele) > 3:
			cmd.load_cgo(invele, self.display_flags['system_name'] + "_" + str(c.idnum) + "_inverted", self.num_frames)
			cmd.load_cgo(numtxt, self.display_flags['system_name'] + "_" + str(c.idnum) + "_invertedindex", self.num_frames)
		bin += 1

  def load_trajectory(self, trajectory_out_fname):

	tbegin = time.time()
	
	#
	# All blobs already have the first frame. They will keep this permanently.
	# All subsequent frames will be readed, loaded, drawn and deleted until failure
	#	

	# Load header and skip first frame (we already have it from the node files)
	try:
		self.traj = FFEA_trajectory.FFEA_trajectory(trajectory_out_fname, load_all = 0, onlyNodes=True)
		try:
			failure = self.traj.skip_frame()
		except:
			failure = 1
	except(IOError):
		failure = 1	

	# Get smallest edge in system
	lmin = float("inf")
	for b in self.blob_list:
		for f in b[0].surf.face:
			l = 2 * f.calc_area(b[0].frames[0])**0.5
			if l < lmin:
				lmin = l

	# Draw first frame
	self.num_frames = 1
	self.draw_frame(self.num_frames - 1, scale = lmin / 20.0)

	# If necessary, stop now (broken traj or user asked for)
	if failure == 1 or self.display_flags['load_trajectory'] != "Trajectory" or self.traj.num_blobs == 0:		
		if failure == 1: 
			print "Failed to load the trajectory: ", failure
		self.wontLoadTraj = 1
		return

	# Else, load rest of trajectory 1 frame at a time, drawing and deleting as we go
	# Save final two frames for later calculations though
	
	while True:
		
		# Get frame from traj
		if self.traj.load_frame(onlyNodes=True) == 0:

			# Scale traj frame
			self.traj.rescale(self.global_scale, -1)
			
			# Load into blob objects asnd increment frame count
			self.add_frame_to_blobs(self.traj)
			self.num_frames += 1

			# Draw whole frame (if above worked, these should work no problem...)
			self.draw_frame(self.num_frames - 1, scale = lmin, draw_static = False)

			# Delete frames from memory
			if(self.num_frames > 3):
				self.traj.delete_frame(index = -3)
	
			self.remove_frame_from_blobs()
		else:
			break

	# Finally show the "progress bar":
	if self.num_frames > 1:
		cmd.mset("1-"+str(self.num_frames))
	# If the trajectory was a single frame, then we loaded nothing:
	else: self.wontLoadTraj = 1

	print "Trajectory loaded in: ", time.time() - tbegin, "s."


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

	# Empty traj object
	self.traj = None

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

	self.waitForTrajToLoad = False # who is closing tkinter window
	
	# Change to any file of names you like
	# fname = os.path.dirname(os.path.realpath(__file__)) + "/system_names_dbzcharacters.txt"
	fname = os.path.dirname(os.path.realpath(__file__)) + "/system_names_greekletters.txt"
	with open(fname, "r") as f:
		for line in f:
			self.system_names.append(line.strip())

	self.display_flags = {
		'matparam': "Plain Solid",
		'show_mesh': "No Mesh",
		'show_numbers': "No Indices", ## PYMOL OK
		'show_pinned': 1,
		'show_beads': "No Beads",
		'show_danger': 0,
		'show_inverted': 1,
		'show_vdw': 0,
		'show_shortest_edge': 0,
		'show_springs': 1,
		'show_box': "No Box",
		'load_trajectory': "Trajectory", ## PYMOL OK
		'highlight': '',
		'load_sfa': 'None',
      'system_name': self.system_names[rint(0, len(self.system_names) - 1)],
      'sele_name': "sele",
      'pin_fname': "",
      'mat_fname': "",
      'which_vdw_type': 0,
      'vdw_type0': -1,
      'vdw_type1': -1,
      'vdw_fname': '',
      'lj_fname': '',
      'lj_eps': 1e15,
      'lj_r0': 1e-9,
      'mat_d': 1.5e3,
      'mat_sm': 370.37e6,
      'mat_bm': 111.11e7,
      'mat_sv': 1e-3,
      'mat_bv': 1e-3}

	self.selected_index = 0
	self.selected_blob = 0
	self.selected_conformation = 0

	self.offset_x = 0
	self.offset_y = 0
	self.offset_z = 0
  
	self.wontLoadTraj = 0 # if traj was not found or there was an error, we'll remember

	# Assume box exists
	self.box_exists = True
	self.box = np.array([-1.0,-1.0,-1.0])
	self.springs = None
	self.lj = None

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

  def draw_frame(self, index, scale = 1.0, draw_static = True):

	# Blobs should only ever have at most 2 frames on them, the initial one and the currently loaded one. So...
	frame_real_index = index

	if index > 0:
		frame_stored_index = 1
	else:
		frame_stored_index = 0
		
	# World first
	if self.display_flags['show_box'] != "No Box":
		if self.box_exists == True:
			self.draw_box(frame_real_index)
		else:
			print "Box does not exist"

	if self.display_flags['show_springs'] == 1 and self.springs != None:
		self.draw_springs(frame_real_index)

	for i in range(self.script.params.num_blobs):
		for j in range(self.script.params.num_conformations[i]):
			frame_real_index = index
			if self.blob_list[i][j].motion_state == "STATIC": 
				if draw_static == True: frame_real_index = "ALL"
				else: continue
			self.blob_list[i][j].draw_frame(frame_stored_index, frame_real_index, self.display_flags, scale = scale)

  def draw_box(self, f):
	
	# A cube has 8 vertices and 12 sides. A hypercube has 16 and 32! "Whoa, that's well cool Ben!" Yeah, ikr 
	obj = [BEGIN, LINES]
	
	# If only outline, no need to loop over entire plane
	#step = self.box
	#step = [b for b in self.box]
	if self.display_flags['show_box'] == "Simulation Box (outline)":
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
			
	elif self.display_flags['show_box'] == "Simulation Box (whole)":
		for i in range(3):
			step = [self.box[i] / self.script.params.es_N[i] for i in range(3)]

		# Loop over the three planes
		for i in range(self.script.params.es_N[0] + 1):
			for j in range(self.script.params.es_N[1] + 1):

				# Get a pair of vertices
				verts = [[i * step[0], j * step[1], 0.0], [i * step[0], j * step[1], self.box[2]]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		# Loop over the three planes
		for i in range(self.script.params.es_N[1] + 1):
			for j in range(self.script.params.es_N[2] + 1):

				# Get a pair of vertices
				verts = [[0.0, i * step[1], j * step[2]], [self.box[0], i * step[1], j * step[2]]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

		# Loop over the three planes
		for i in range(self.script.params.es_N[0] + 1):
			for j in range(self.script.params.es_N[2] + 1):

				# Get a pair of vertices
				verts = [[i * step[0], 0.0, j * step[2]], [i * step[0], self.box[1], j * step[2]]]
				
				for l in range(2):
					obj.extend([VERTEX, verts[l][0], verts[l][1], verts[l][2]])

				
					

	obj.append(END)
	cmd.load_cgo(obj, self.display_flags['system_name'] +"_Simulation_Box", f + 1)

  def draw_springs(self, f):

      for s in self.springs.spring:
	# print(self.springs.spring.index(s))
         # Get correct frames
         correct_frame = [-1 for i in range(self.script.params.num_blobs)]
         for i in range(self.script.params.num_blobs):
            if self.blob_list[i][0].motion_state == "STATIC":
               correct_frame[i] = 0
         # print "correct_frame: ", correct_frame

         # Draw, because this spring exists
         try:
	  # s.print_details()
           springjoints = np.array([self.blob_list[s.blob_index[i]][s.conformation_index[i]].frames[correct_frame[s.blob_index[i]]].pos[s.node_index[i]][0:3] for i in range(2)])
	   #print(springjoints)
         except(AttributeError):
         #  print("Whut!")
           continue
         except:
           print "Something went wrong with this spring"
		#	except(IndexError):
			#	if s.blob_index[i] >= self.num_blobs:
			#		print "fuck"
			#	if s.conformation_index[i] >= self.num_conformations[i]:
			#		print "fuck2"
			#	if s.node_index[i] >= self.blob_list[i].

         # Axes for helix
         zax = springjoints[1] - springjoints[0]
         l = np.linalg.norm(zax)
         zax = zax / l

	 if np.fabs(np.dot(zax, np.array([1.0,0.0,0.0]))) < 1.0:
                 xax = np.cross(zax, np.array([1.0,0.0,0.0]))
                 yax = np.cross(zax, xax)
	 else:
	         xax = np.cross(zax, np.array([0.0,1.0,0.0]))
	         yax = np.cross(zax, xax)

         xax = xax / np.linalg.norm(xax)
         yax = yax / np.linalg.norm(yax)

         # Radius of helix (let original radius be 5A, poisson ratio = 0.01)
         r = 2 - 0.01 * (l - s.l)

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
         


