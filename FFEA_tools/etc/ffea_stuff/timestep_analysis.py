#!/usr/bin/env python

import os, sys, subprocess
from math import *

def calc_decimal_places(a):
	i = 0
	while a != 0:
		a = a - int(a)
		a = a * 10
		i = i + 1
		print a
	return i - 1
		

def set_scale(scale, inputffea):
	infile = open(inputffea, "r")
	lines = infile.readlines()
	infile.close()
	outfile = open(inputffea, "w")
	for line in lines:
		if line == "\n":
			outfile.write("\n")
		elif line.split()[0] == "<scale":
			outfile.write("\t\t<scale = " + str(scale) + ">\n")	
		else:
			outfile.write(line)
	outfile.close()

def set_dt(dt, inputffea):
	infile = open(inputffea, "r")
	lines = infile.readlines()
	infile.close()
	outfile = open(inputffea, "w")
	for line in lines:
		if line == "\n":
			outfile.write("\n")
		elif line.split()[0] == "<dt":
			outfile.write("\t<dt = " + str(dt) + ">\n")	
		else:
			outfile.write(line)	
	outfile.close()

if len(sys.argv) != 5:
	sys.exit("Usage: python timestep_analysis.py [INPUT .FFEA FILE] [MIN SCALE] [MAX SCALE] [NUM DATA POINTS]")

input_ffea = sys.argv[1]
min_scale = float(sys.argv[2])
max_scale = float(sys.argv[3])
num_data_points = int(sys.argv[4])

output_data = open("timestep_stuff.csv", "w")
output_data.write("Scale,Max_timestep\n\n")
scales = [0.0 for i in range(num_data_points)]


for i in range(len(scales)):
	scales[i] = min_scale * pow(max_scale/min_scale, i/float(num_data_points))

print scales
for scale in scales:
	dt_upper = 1e-8
	dt_lower = 1e-15
	dt = 1e-8
	set_scale(scale, input_ffea)
	set_dt(dt, input_ffea)
	iterations = 0
	while iterations < 70:
		count = 0
		while count < 1:
			if(subprocess.call(["./ffea",input_ffea]) == 0):
				count = count + 1
			else:
				break
		if count == 1:
			dt_lower = dt
		else:
			dt_upper = dt
		
		dt = 0.5 * (dt_lower + dt_upper)
		set_dt(dt, input_ffea)
		iterations = iterations + 1
	output_data.write(str(scale) + "," + str(dt) + "\n")
				
