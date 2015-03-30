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

import FFEA_viewer_control_window
import FFEA_viewer_display_window


file_to_load = None
energy_thresh = 1e-6
num_frames_to_read = float("inf")

for i in range(1, len(sys.argv), 2):
	if sys.argv[i] == "-script":
		try:
			file_to_load = sys.argv[i + 1]
		except IndexError:
			print "Parameter not specified"
			break

	elif sys.argv[i] == "-frames":
		try:
			num_frames_to_read = int(sys.argv[i + 1])
		except IndexError:
			print "Parameter not specified"
			break
	elif sys.argv[i] == "-energy":
		try:
			energy_thresh = float(sys.argv[i + 1])
		except IndexError:
			print "Parameter not specified"
			break
	else:
		print "Unrecognised flag '" + sys.argv[i] + "'\n"
		break


root = Tk()
root.geometry("650x600+801+30")
root.title("FFEA Viewer - Control")
control_window = FFEA_viewer_control_window.FFEA_viewer_control_window(root, file_to_load, num_frames_to_read, energy_thresh=energy_thresh)
root.mainloop()
