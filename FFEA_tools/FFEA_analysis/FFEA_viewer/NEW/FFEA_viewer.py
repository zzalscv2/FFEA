from Tkinter import *

import tkFileDialog
import tkColorChooser

from multiprocessing import Process, Pipe

import sys
import os
import random
import math
import time

import FFEA_viewer_control_window
#import FFEA_viewer_display_window

def print_help():
	print "Usage: python " + sys.argv[0] + " [options]"
	print "Options:"
	print "  -h/--help\t\t\t\tdisplay this help"
	print "  -s/--script 'fname'\t\t\tscript fname to use"
	print "  -e/--energy 'energy_thresh'\t\tenergy threshold (no idea what this does!)"

file_to_load = None
energy_thresh = 1e-6
num_frames_to_read = float("inf")

for i in range(1, len(sys.argv), 2):
	if os.path.exists(sys.argv[i]) or os.path.exists(os.path.abspath(sys.argv[i])):
		file_to_load = sys.argv[i]
		break

	elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
		print_help()
		sys.exit()
	elif sys.argv[i] == "-s" or sys.argv[i] == "--script":
		try:
			file_to_load = sys.argv[i + 1]
		except IndexError:
			print "Parameter not specified"
			continue

	elif sys.argv[i] == "-e" or sys.argv[i] == "--energy":
		try:
			energy_thresh = float(sys.argv[i + 1])
		except IndexError:
			print "Parameter not specified"
			continue
	else:
		print "Unrecognised flag '" + sys.argv[i] + "'\n"
		break


root = Tk()
root.geometry("650x600+801+30")
root.title("FFEA Viewer - Control")
control_window = FFEA_viewer_control_window.FFEA_viewer_control_window(root, file_to_load, energy_thresh=energy_thresh)
root.mainloop()
