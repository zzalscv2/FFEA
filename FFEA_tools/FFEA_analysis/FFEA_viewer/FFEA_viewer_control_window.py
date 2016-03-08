from Tkinter import *

import tkFileDialog
import tkColorChooser

from multiprocessing import Process, Pipe

import sys
import os
import random
import math
import time

#import pyximport; pyximport.install()
import FFEA_viewer_display_window

class FFEA_viewer_control_window:

	def __init__(self, root, file_to_load, num_frames_to_read, energy_thresh=1.0e6):
		self.energy_threshold = energy_thresh
		self.num_frames_to_read = num_frames_to_read
		self.master = root
		self.master.protocol("WM_DELETE_WINDOW", self.death)
		self.master.bind('<Escape>', self.esc_key_handler)
		
		frame = Frame(self.master)
		frame.pack()
		
		menubar = Menu(frame)
		
		filemenu = Menu(menubar, tearoff=0)
		filemenu.add_command(label="Load 'ffea' file", command=self.choose_ffea_file_to_load)
		filemenu.add_command(label="Import...", command=self.import_handler)
		
		emdbmenu = Menu(menubar, tearoff=0)
		emdbmenu.add_command(label="EMDB", command=self.emdb_wizard)
		
		menubar.add_cascade(label="File", menu=filemenu)
		menubar.add_cascade(label="EMDB", menu=emdbmenu)
		
		self.master.config(menu=menubar)
		
		# Top frame
		top_frame = Frame(self.master)
		top_frame.pack()
		
		# scrollable list box of loaded blobs
		listbox_frame = Frame(top_frame)
		listbox_frame.pack(side=LEFT, anchor=N)
		listbox_scrollbar = Scrollbar(listbox_frame)
		listbox_scrollbar.pack(side=RIGHT, fill=Y)
		self.blob_listbox = Listbox(listbox_frame, width=40, height=10, yscrollcommand=listbox_scrollbar.set, selectmode=SINGLE)
		listbox_scrollbar.config(command=self.blob_listbox.yview)
		self.blob_listbox.pack(fill=Y)
		self.blob_listbox.bind('<<ListboxSelect>>', self.on_blob_listbox_select)
		
		# Selected blob frame
		selected_index_frame = Frame(top_frame, relief=RAISED, bd=1, bg='#ffffc2')
		selected_index_frame.pack(side=LEFT, fill=Y)
		self.selected_index_label = Label(selected_index_frame, text="No Blob selected...", justify=LEFT, bg='#ffffc2')
		self.selected_index_label.pack()

		self.show_hide = IntVar()
		self.check_button_show_hide = Checkbutton(selected_index_frame, text="Show/Hide Blob", variable=self.show_hide, command=self.hide_blob)
		self.check_button_show_hide.pack(side=TOP, anchor=N)
		vdw_frame = Frame(selected_index_frame, bd=1, bg='#ffffc2')
		vdw_frame.pack(side=TOP, fill=X)
		binding_frame = Frame(selected_index_frame, bd=1, bg='#ffffc2')
		binding_frame.pack(side=TOP, fill=X)

		self.edit_vdw = IntVar()
		self.check_button_edit_vdw = Checkbutton(vdw_frame, text="Edit VdW", variable=self.edit_vdw, state=DISABLED, command=self.something_has_changed)
		self.check_button_edit_vdw.pack(side=LEFT)
		self.save_button_vdw = Button(vdw_frame, text="Save VdW...", state=DISABLED, command=self.save_vdw)
		self.save_button_vdw.pack(side=RIGHT)
		self.edit_binding_sites = IntVar()
		self.check_button_edit_binding = Checkbutton(binding_frame, text="Edit Binding Site:", variable=self.edit_binding_sites, state=DISABLED, command=self.something_has_changed)
		self.check_button_edit_binding.pack(side=LEFT)
		self.save_button_binding = Button(binding_frame, text="Save Binding Sites...", state=DISABLED)#, command=self.save_binding_sites)
		self.save_button_binding.pack(side=RIGHT)
		self.bsite_spin = Spinbox(binding_frame, state=DISABLED)
		self.bsite_spin.pack(side=LEFT)
		self.blob_info_list = []
		
		# Display flags frame
		display_flags_frame = Frame(self.master, relief=SUNKEN, bd=1)
		display_flags_frame.pack(side=RIGHT, anchor=N)
		self.show_mesh = IntVar()
		check_button_show_mesh = Checkbutton(display_flags_frame, text="Mesh", variable=self.show_mesh, command=self.something_has_changed)
		check_button_show_mesh.pack(anchor=W)
		self.show_mesh_surf = IntVar()
		check_button_show_mesh_surf = Checkbutton(display_flags_frame, text="Mesh surf", variable=self.show_mesh_surf, command=self.something_has_changed)
		check_button_show_mesh_surf.pack(anchor=W)
		self.show_solid = IntVar(value=1)
		check_button_show_solid = Checkbutton(display_flags_frame, text="Solid", variable=self.show_solid, command=self.something_has_changed)
		check_button_show_solid.pack(anchor=W)
		self.show_flat = IntVar(value=0)
		check_button_show_flat = Checkbutton(display_flags_frame, text="Flat", variable=self.show_flat, command=self.something_has_changed)
		check_button_show_flat.pack(anchor=W)
		self.show_material = IntVar(value=0)
		check_button_show_material = Checkbutton(display_flags_frame, text="Show Material", variable=self.show_material, command=self.something_has_changed)
		check_button_show_material.pack(anchor=W)
		self.show_vdw_only = IntVar(value=0)
		check_button_show_vdw_only = Checkbutton(display_flags_frame, text="VdW only", variable=self.show_vdw_only, command=self.something_has_changed)
		check_button_show_vdw_only.pack(anchor=W)
		self.show_shortest_edge = IntVar(value=0)
		check_button_show_shortest_edge = Checkbutton(display_flags_frame, text="Shortest edge", variable=self.show_shortest_edge, command=self.something_has_changed)
		check_button_show_shortest_edge.pack(anchor=W)
		self.show_inverted = IntVar(value=0)
		check_button_show_inverted = Checkbutton(display_flags_frame, text="Inverted elements", variable=self.show_inverted, command=self.something_has_changed)
		check_button_show_inverted.pack(anchor=W)
		self.show_node_numbers = IntVar()
		check_button_show_node_numbers = Checkbutton(display_flags_frame, text="Node numbers", variable=self.show_node_numbers, command=self.something_has_changed)
		check_button_show_node_numbers.pack(anchor=W)
		self.show_linear_nodes = IntVar()
		check_button_show_linear_nodes = Checkbutton(display_flags_frame, text="Linear only", variable=self.show_linear_nodes, command=self.something_has_changed)
		check_button_show_linear_nodes.pack(anchor=W)
		self.show_box = IntVar(value=0)
		check_button_show_box = Checkbutton(display_flags_frame, text="Box", variable=self.show_box, command=self.something_has_changed)
		check_button_show_box.pack(anchor=W)
		self.hide_frozen = IntVar()
		check_button_hide_frozen = Checkbutton(display_flags_frame, text="Hide Frozen", variable=self.hide_frozen, command=self.something_has_changed)
		check_button_hide_frozen.pack(anchor=W)
		self.blob_colour = (0.0,0.588,1.0)
		self.choose_colour = Button(display_flags_frame, text="Colour", command=self.choose_blob_colour, bg=self.get_colour_code(self.blob_colour))
		self.choose_colour.pack(anchor=W)
		
		self.recording = IntVar()
		check_button_recording = Checkbutton(display_flags_frame, text="Recording", variable=self.recording, command=self.something_has_changed)
		check_button_recording.pack(anchor=W, )
		
		# Pause loading button
		self.pause_loading_button = Button(self.master, text="Pause loading", command=self.pause_loading_handler, state=DISABLED)
		self.pause_loading_button.pack()
		self.pause_loading = False

		# Frame slider
		frame_slider_frame = Frame(self.master)
		frame_slider_frame.pack()
		self.slider_current_frame = IntVar(value=0)
		self.frame_number_label = Label(frame_slider_frame, text="0/0")
		self.frame_number_label.pack()
		frame_label = Label(frame_slider_frame, text="Frame")
		frame_label.pack(side=LEFT)
		self.frame_slider = Scale(frame_slider_frame, orient=HORIZONTAL, variable=self.slider_current_frame, length=250, to=0, width=1, showvalue=0, command=self.frame_slide_handler)
		self.frame_slider.pack(side=LEFT)
		
		# Start/stop button
		self.start_stop_button = Button(frame_slider_frame, text=">", command=self.start_stop_handler, state=DISABLED)
		self.start_stop_button.pack(side=RIGHT)
		
		# Speed slider
		speed_slider_frame = Frame(self.master)
		speed_slider_frame.pack()
		self.slider_speed = IntVar(value=1)
		self.speed_label = Label(speed_slider_frame, text="Speed")
		self.speed_label.pack(side=LEFT)
		self.speed_slider = Scale(speed_slider_frame, orient=HORIZONTAL, variable=self.slider_speed, length=100, from_=1,to=10, width=1, showvalue=0, command=self.speed_slide_handler)
		self.speed_slider.pack(side=RIGHT)
		
		# Projection type frame
		projection_frame = Frame(self.master, relief=SUNKEN, bd=1)
		projection_frame.pack()
		self.projection = StringVar()
		rb_persp = Radiobutton(projection_frame, text="Perspective", variable=self.projection, value="perspective", command=self.something_has_changed)
		rb_ortho = Radiobutton(projection_frame, text="Orthographic", variable=self.projection, value="orthographic", command=self.something_has_changed)
		rb_persp.pack(side=LEFT)
		rb_ortho.pack(side=LEFT)
		rb_ortho.select()

		# screenshot button
		self.screenshot_button = Button(self.master, text="Screenshot", command=self.screenshot)
		self.screenshot_button.pack()

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
		self.master.after(100, self.get_updates_from_display())

		# If a file name was given, open it
		if file_to_load != None:
			self.load_ffea(file_to_load)
		
	def hide_blob(self):
		if self.show_hide.get() == 1:
			self.speak_to_display.send({'hide_blob':self.selected_index})
		else:
			self.speak_to_display.send({'show_blob':self.selected_index})

	def screenshot(self):
		# set up the options for the save file dialog box
		options = {}
		options['defaultextension'] = '.tga'
		options['filetypes'] = [('TGA files', '.tga'), ('all files', '.*')]
		options['initialdir'] = os.getcwd()
		options['title'] = 'Save TGA file'

		# Ask user to select a file
		chosen_tga_name = tkFileDialog.asksaveasfilename(**options)
		if len(chosen_tga_name) == 0:
			return

		self.speak_to_display.send({'save_screenshot':chosen_tga_name})

	def save_vdw(self):
		# set up the options for the save file dialog box
		options = {}
		options['defaultextension'] = '.vdw'
		options['filetypes'] = [('vdw files', '.vdw'), ('all files', '.*')]
		options['initialdir'] = os.getcwd()
		options['title'] = 'Save vdw file'

		# Ask user to select a file
		chosen_vdw_name = tkFileDialog.asksaveasfilename(**options)
		if len(chosen_vdw_name) == 0:
			return

		self.speak_to_display.send({'save_vdw':chosen_vdw_name})

	def choose_blob_colour(self):
		new_colour = tkColorChooser.askcolor(initialcolor=self.get_colour_code(self.blob_colour), title ="Choose Blob display colour")
		if new_colour[0] == None:
			return
		new_colour = new_colour[0]
		print "Chose new blob colour:", new_colour
		self.blob_colour = (new_colour[0]/255.0, new_colour[1]/255.0, new_colour[2]/255.0)
		self.choose_colour.config(bg=self.get_colour_code(self.blob_colour))
		self.something_has_changed()

	def on_blob_listbox_select(self, event):
		i = int(event.widget.curselection()[0])
		blob_info_text = "Selected Blob " + str(i) + "\n"
		blob_info_text += self.blob_info_list[i]
		self.selected_index_label.config(text=blob_info_text)
		self.save_button_vdw.config(state=NORMAL)
		self.check_button_edit_vdw.config(state=NORMAL)
		self.save_button_binding.config(state=NORMAL)
		self.check_button_edit_binding.config(state=NORMAL)
		self.change_indices(i)
		self.something_has_changed()

	def change_indices(self, index):

		# Change to control blob index to selected_blob and selected_conformation
		self.selected_index = index
		k = 0
		for i in range(self.num_blobs):
			for j in range(self.num_conformations[i]):
				if self.selected_index == k:
					self.selected_blob = i
					self.selected_conformation = j
					return
				else:
					k += 1

	def launch_display_window(self, speak_to_control, ffea_fname):
		FFEA_viewer_display_window.FFEA_viewer_display_window(speak_to_control, ffea_fname, self.num_frames_to_read, energy_thresh=self.energy_threshold)

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

	def load_ffea(self, ffea_fname):

		# Kill any existing display window and reset
		if self.display_window_exists == True:
			self.send_message_to_display(True)
			time.sleep(1)
			self.reset_all()

		# Check if given file exists
		if os.path.isfile(ffea_fname) == False:
			print "No such file:", ffea_fname
			return

		# Launch the display window, with a pipe set up so the two processes can "speak" to each other
		self.speak_to_display, speak_to_control = Pipe()
		self.display_window_process = Process(target=launch_display_window, args=(speak_to_control, ffea_fname, self.num_frames_to_read, self.energy_threshold))
		self.display_window_process.start()
		self.display_window_exists = True

		self.start_stop_button.config(text=">")
		self.pause_loading_button.config(text="Pause loading")
		
		self.frame_slider.config(width=15)
		self.start_stop_button.config(state=NORMAL)
		self.pause_loading_button.config(state=NORMAL)
		self.speed_slider.config(width=10)
		
		self.there_is_something_to_send_to_display_window = True
		
	def send_message_to_display(self, death):
		self.speak_to_display.send({	'animate':self.animate,
						'pause_loading':self.pause_loading,
						'display_flags':{	'show_mesh': self.show_mesh.get(),
									'show_solid': self.show_solid.get(),
									'show_flat': self.show_flat.get(),
									'show_material': self.show_material.get(),
									'show_vdw_only': self.show_vdw_only.get(),
									'show_node_numbers': self.show_node_numbers.get(),
									'vdw_edit_mode': self.edit_vdw.get(),
									'binding_site_edit_mode': self.edit_binding_sites.get(),
									'blob_colour': self.blob_colour,
									'selected_index': self.selected_index,
									'selected_blob': self.selected_blob,
									'selected_conformation': self.selected_conformation,
									'show_shortest_edge': self.show_shortest_edge.get(),
									'show_inverted': self.show_inverted.get(),
									'show_pinned_nodes': 1,
									'show_linear_nodes_only': self.show_linear_nodes.get(),
									'hide_frozen': self.hide_frozen.get(),
									'show_mesh_surf': self.show_mesh_surf.get()
								},
						'speed': self.slider_speed.get(),
						'show_box': self.show_box.get(),
						'change_frame_to': self.change_frame_to,
						'selected_index': self.selected_index,
						'recording': self.recording.get(),
						'projection': self.projection.get(),
						'death': death})

	def get_colour_code(self, rgb_float_tuple):
		c = (rgb_float_tuple[0] * 255, rgb_float_tuple[1] * 255, rgb_float_tuple[2] * 255)
		return '#%02x%02x%02x' % c

	def import_handler(self):
		print "Import"

	def something_has_changed(self):
		self.there_is_something_to_send_to_display_window = True

	def frame_slide_handler(self, f):
		self.change_frame_to = f
		self.something_has_changed()

	def speed_slide_handler(self, s):
		self.something_has_changed()

	def emdb_wizard(self):
		print "yo"
	
	def start_stop_handler(self):
		if self.animate == True:
			self.animate = False
			self.start_stop_button.config(text=">")
		else:
			self.animate = True
			self.start_stop_button.config(text="X")
		self.something_has_changed()

	def pause_loading_handler(self):
		if self.pause_loading == True:
			self.pause_loading = False
			self.pause_loading_button.config(text="Pause loading")
		else:
			self.pause_loading_button.config(text="Pausing...", state=DISABLED)
			self.pause_loading = True
		self.something_has_changed()

	def esc_key_handler(self, event):
		print "Escape key pressed."
		self.death()

	def death(self):
		if self.display_window_exists == True:
			print "Telling display to die..."
			self.send_message_to_display(True)
		sys.exit("Control Window: Adios.")
	
	def get_updates_from_display(self):
		if self.display_window_exists == True:
			if self.speak_to_display.poll() == True:
				display_stuff = self.speak_to_display.recv()

				if 'num_blobs' in display_stuff.keys():
					self.num_blobs = display_stuff['num_blobs']
				if 'num_conformations' in display_stuff.keys():
					self.num_conformations = display_stuff['num_conformations']

				if "add_blob" in display_stuff.keys():
					add_blob_info = display_stuff['add_blob']
					self.blob_listbox.insert(END, add_blob_info['name'])
					self.blob_info_list.append(add_blob_info['info'])
					self.master.after(100, self.get_updates_from_display())
					return

				if display_stuff['death'] == True:
					print "Display window closed."
					self.reset_all()
					self.master.after(100, self.get_updates_from_display())
					return

				self.frame_slider.config(to=display_stuff['num_frames'])
				if display_stuff['current_frame'] != -1:
					self.slider_current_frame.set(display_stuff['current_frame'])
					self.frame_number_label.config(text=str(display_stuff['current_frame']) + "/" + str(display_stuff['num_frames']))
				else:
					self.frame_number_label.config(text=str(self.slider_current_frame.get()) + "/" + str(display_stuff['num_frames']))

				if self.pause_loading == True and display_stuff['pausing'] == False:
					self.pause_loading_button.config(text="Pausing...", state=DISABLED)
				else:
					if self.pause_loading == True:
						self.pause_loading_button.config(text="Resume loading", state=NORMAL)
					else:
						self.pause_loading_button.config(text="Pause loading", state=NORMAL)

			if self.there_is_something_to_send_to_display_window == True:
				self.there_is_something_to_send_to_display_window = False
				self.send_message_to_display(False)
				self.change_frame_to = -1

		self.master.after(100, self.get_updates_from_display)

	def reset_all(self):
		self.display_window_exists = False
		self.animate = False
		self.there_is_something_to_send_to_display_window = False
		self.selected_index = 0
		self.frame_slider.config(to=0, width=1)
		self.speed_slider.config(width=1)
		self.slider_current_frame.set(0)
		self.frame_number_label.config(text="0/0")

		self.blob_listbox.delete(0,self.blob_listbox.size())

		self.pause_loading = False
		self.pause_loading_button.config(text="Pause loading", state=DISABLED)

		self.blob_info_list = []
		self.selected_index_label.config(text="No Blob selected...")
		self.save_button_vdw.config(state=DISABLED)
		self.check_button_edit_vdw.config(state=DISABLED)
		self.save_button_binding.config(state=DISABLED)
		self.check_button_edit_binding.config(state=DISABLED)
		self.edit_vdw.set(0)
		self.edit_binding_sites.set(0)
		self.speak_to_display.close()

def launch_display_window(speak_to_control, ffea_fname, num_frames_to_read, energy_threshold):
	FFEA_viewer_display_window.FFEA_viewer_display_window(speak_to_control, ffea_fname, num_frames_to_read, energy_thresh=energy_threshold)