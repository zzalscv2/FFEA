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
if len(sys.argv) > 1:
	file_to_load = sys.argv[1]


if len(sys.argv) == 3:
	energy_thresh = float(sys.argv[2])
else:
	energy_thresh = 1.0e6


root = Tk()
root.geometry("650x600+801+30")
root.title("FFEA Viewer - Control")
control_window = FFEA_viewer_control_window.FFEA_viewer_control_window(root, file_to_load, energy_thresh=energy_thresh)
root.mainloop()
