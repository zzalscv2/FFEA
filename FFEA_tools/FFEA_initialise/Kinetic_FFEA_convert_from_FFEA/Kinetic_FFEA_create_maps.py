#!/usr/bin/env python

import sys, os
from math import *
from Vectors import vector3
import FFEA_traj

def get_nodes_from_script(fname):

	# Get fname content
	fin = open(fname, "r")
	lines = fin.readlines()
	fin.close()
	
	# Get node fnames (should be 2)
	node_fnames = []
	for line in lines:
	
		# Break if no equals
		if "=" not in line:
			continue
		
		line = line.strip()[1:-1]
		lvalue, rvalue = line.split("=")
		if lvalue.strip() == "nodes":
			node_fnames.append(rvalue.strip())
			
	if len(node_fnames) != 2:
		sys.exit("Error. In " + fname + " there should be a two blobs, each with a single conformation giving two node filenames. You have " + str(len(node_fnames)) + "filenames\n")
	else:
		return node_fnames

def get_topologies_from_script(fname):

	# Get fname content
	fin = open(fname, "r")
	lines = fin.readlines()
	fin.close()
	
	# Get top fnames (should be 2)
	top_fnames = []
	for line in lines:
	
		# Break if no equals
		if "=" not in line:
			continue
		
		line = line.strip()[1:-1]
		lvalue, rvalue = line.split("=")
		if lvalue.strip() == "topology":
			top_fnames.append(rvalue.strip())
			
	if len(top_fnames) != 2:
		sys.exit("Error. In " + fname + " there should be a two blobs, each with a single conformation giving two topology filenames. You have " + str(len(top_fnames)) + "filenames\n")
	else:
		return top_fnames

def get_traj_from_script(fname):

	# Get fname content
	fin = open(fname, "r")
	lines = fin.readlines()
	fin.close()
	
	# Get traj fname
	node_fnames = []
	for line in lines:
	
		# Break if no equals
		if "=" not in line:
			continue
		
		line = line.strip()[1:-1]
		lvalue, rvalue = line.split("=")
		if lvalue.strip() == "trajectory_out_fname":
			return rvalue.strip()

def get_scale_from_script(fname):

	# Get fname content
	fin = open(fname, "r")
	lines = fin.readlines()
	fin.close()
	
	# Get traj fname
	node_fnames = []
	for line in lines:
	
		# Break if no equals
		if "=" not in line:
			continue
		
		line = line.strip()[1:-1]
		lvalue, rvalue = line.split("=")
		if lvalue.strip() == "scale":
			return float(rvalue.strip())

def get_num_steps_from_script(fname):

	# Get fname content
	fin = open(fname, "r")
	lines = fin.readlines()
	fin.close()
	
	# Get traj fname
	node_fnames = []
	for line in lines:
	
		# Break if no equals
		if "=" not in line:
			continue
		
		line = line.strip()[1:-1]
		lvalue, rvalue = line.split("=")
		if lvalue.strip() == "num_steps":
			return float(rvalue.strip())

def make_a_spring_file(fname, springs, k, l):

	fout = open(fname, "w")
	fout.write("ffea spring file\nnum_springs " + str(len(springs)) + "\nblob0 conformation0 node0 blob1 conformation1 node1 k l:\n")
	for spring in springs:
		fout.write("%d %d %d %d %d %d %e %e\n" % (0, 0, spring[0], 1, 0, spring[1], k, l))

	fout.close()

def make_initial_script(orig_script, new_script, orig_node_fnames, new_node_fname, new_spring_fname, spring_k, spring_l):

	fin = open(orig_script, "r")
	fout = open(new_script, "w")
	lines = fin.readlines()

	num_blobs = 0
	one_static = 0
	node_fnames_read = 0
	for inline in lines:
		if inline.strip()[1:-1] == "spring" or inline.strip()[1:-1] == "/spring":
			continue

		if inline.strip()[1:-1] == "/system":
			fout.write("\t<spring>\n")
			fout.write("\t\t<spring_fname = %s>\n" % (spring_fname))
			fout.write("\t</spring>\n")

		elif "=" not in inline:
			fout.write(inline)
			continue
		
		sline = inline.strip()[1:-1].split("=")
		outline = inline
		if sline[0].strip() == "spring_fname":
			continue
		elif sline[0].strip() == "calc_noise":
			outline = "\t<calc_noise = 0>\n"
		elif sline[0].strip() == "num_blobs":
			if sline[1].strip() != "2":
				sys.exit("Error. Specified " + sline[1].strip() + " blobs. There must be two to create a mapping.\n")
		elif sline[0].strip() == "restart":
			outline = "\t<restart = 0>\n"
		elif sline[0].strip() == "num_steps":
			outline = "\t<num_steps = 1e3>\n"
		elif sline[0].strip() == "check":
			outline = "\t<check = 1e2>\n"
		elif sline[0].strip() == "num_conformations":
			outline = "\t<num_conformations = (1,1)>\n"
		elif sline[0].strip() == "num_states":
			outline = "\t<num_states = (1,1)>\n"
		elif sline[0].strip() == "motion_state":
			if one_static == 0:
				outline = "\t\t\t<motion_state = DYNAMIC>\n"
			else:
				outline = "\t\t\t<motion_state = DYNAMIC>\n"

			one_static += 1

		elif sline[0].strip() == "nodes":
			node_fnames_read += 1
			if node_fnames_read == 1:
				outline = "\t\t\t<nodes = " + new_node_fname + ">\n"
		
		fout.write(outline)


	fin.close()
	fout.close()

# Change restart in .ffea script
def change_restart(fname):
	
	fin = open(fname, "r")
	fout = open("temp.ffea", "w")
	lines = fin.readlines()
	fin.close()
	for line in lines:
		if "=" not in line:
			fout.write(line)
		else:
			lvalue, rvalue = line.strip()[1:-1].split("=")	
			if lvalue.strip() == "restart":
				fout.write("\t<restart = 1>\n")
			else:
				fout.write(line)
	
	fin.close()
	fout.close()
	os.system("cp temp.ffea " + fname)
	os.system("rm temp.ffea")

# Add another round of simulation on
def add_sim_time(fname, time):
	
	fin = open(fname, "r")
	fout = open("temp.ffea", "w")
	lines = fin.readlines()
	fin.close()
	for line in lines:
		if "=" not in line:
			fout.write(line)
		else:
			lvalue, rvalue = line.strip()[1:-1].split("=")	
			if lvalue.strip() == "num_steps":
				fout.write("\t<num_steps = %e>\n" % (float(rvalue) + time))
			else:
				fout.write(line)
	
	fin.close()
	fout.close()
	os.system("cp temp.ffea " + fname)
	os.system("rm temp.ffea")

# Convert final frame of trajectory to a .node file
def final_frame_to_nodes(traj_fname, moving_node_fname, scale):
	
	# Get stuff from .node file first
	fin = open(moving_node_fname, "r")
	fin.readline()
	num_nodes = int(fin.readline().split()[1])	
	num_surface_nodes = int(fin.readline().split()[1])
	num_interior_nodes = int(fin.readline().split()[1])
	fin.close()

	# Find final frame in .traj file firstly getting the number of frames
	frame_count = 0
	num_frames = FFEA_traj.get_num_frames(traj_fname)
	fin = open(traj_fname, "r")
	while True:
		if fin.readline().strip() == "*":
			frame_count += 1
			if frame_count == num_frames:
				break
	
	# Moving fname corresponds to blob 0
	fin.readline()
	fin.readline()
	fout = open(moving_node_fname, "w")
	fout.write("ffea node file\n")
	fout.write("num_nodes %d\n" % (num_nodes))
	fout.write("num_surface_nodes %d\n" % (num_surface_nodes))
	fout.write("num_interior_nodes %d\n" % (num_interior_nodes))

	# Surface nodes
	fout.write("surface nodes:\n")
	for i in range(num_surface_nodes):
		sline = fin.readline().split()
		fout.write("%6.3f %6.3f %6.3f\n" % (float(sline[0]) / scale, float(sline[1]) / scale, float(sline[2]) / scale))

	# Interior nodes
	fout.write("interior nodes:\n")
	for i in range(num_interior_nodes):
		sline = fin.readline().split()
		fout.write("%6.3f %6.3f %6.3f\n" % (float(sline[0]) / scale, float(sline[1]) / scale, float(sline[2]) / scale))
	
	fin.close()
	fout.close()

def make_nodes_from_traj(final_node_fnames, orig_node_fnames, traj_fname, num_frames):
	
	# Get traj from file
 	traj = FFEA_traj.FFEA_traj(traj_fname, num_frames)
	
	for i in range(len(orig_node_fnames)):
		fin = open(orig_node_fnames[i], "r")
		fout = open(final_node_fnames[i], "w")
		
		fout.write(fin.readline())

		# num_nodes			
		line = fin.readline()
		num_nodes = int(line.split()[1])
		fout.write(line)

		# num_surface_nodes			
		line = fin.readline()
		num_surface_nodes = int(line.split()[1])
		fout.write(line)

		# num_interior_nodes			
		line = fin.readline()
		num_interior_nodes = int(line.split()[1])
		fout.write(line)

		fin.close()

		# surface nodes
		fout.write("surface nodes:\n")
		for j in range(num_surface_nodes):
			fout.write("%6.3e %6.3e %6.e3\n" % (traj.blob[i].frame[-1].node[j].x, traj.blob[i].frame[-1].node[j].y, traj.blob[i].frame[-1].node[j].z))

		# interior nodes
		fout.write("interior nodes:\n")
		for j in range(num_surface_nodes, num_nodes):
			fout.write("%6.3e %6.3e %6.e3\n" % (traj.blob[i].frame[-1].node[j].x, traj.blob[i].frame[-1].node[j].y, traj.blob[i].frame[-1].node[j].z))

# Main program begins here
if len(sys.argv) < 3:
	sys.exit("Usage python " + sys.argv[0] + " [INPUT .ffea script (must have 2 blobs)] [tolerance length] [Pairs of nodes {[i,j]}]\n")
	
# Get initial scripts
script_fname = sys.argv[1]
path = os.path.dirname(script_fname)
if path == '':
	path = '.'
script_fname = path + "/" + script_fname
script_fname_base = script_fname.split("/")[1].split(".")[0]

# Working script
initial_script_fname = script_fname_base + "_initial.ffea"

# Get node fnames from script
node_fname = get_nodes_from_script(script_fname)

# Get top fnames from script
top_fname = get_topologies_from_script(script_fname)

# Get traj fname from script
traj_fname = get_traj_from_script(script_fname)

# Get scale_from_script
scale = get_scale_from_script(script_fname)

# Make working node files
moving_node_fname = script_fname_base + "_working.node"
os.system("cp " + node_fname[0] + " " + moving_node_fname)

# Spring constants for these springs
l = float(sys.argv[2])
k = 1.5e-2

# Get node pairs for springs
node_pairs = []
for i in range(3, len(sys.argv), 1):
	argsplit = sys.argv[i][1:-1].strip().split(",")
	try:
		node_pair = [int(argsplit[0]), int(argsplit[1])]
	except(ValueError):
		sys.exit("Error. node_pair %d not correctly formatted. Format is [i,j]\n" % (i - 3))

	node_pairs.append(node_pair)	
	
num_pairs = len(node_pairs)
spring_fname = script_fname_base + "_mapping.spring"

# Make initial script
make_initial_script(script_fname, initial_script_fname, node_fname, moving_node_fname, spring_fname, k, l)
make_a_spring_file(spring_fname, node_pairs, k, l)

# Get num_steps_from_script
num_steps = get_num_steps_from_script(initial_script_fname)

# Enter minimisation loop
not_satisfied = True
num_runs = 0
while not_satisfied:
	
	if num_runs == 0:
		os.system("python ../Node_tools/translate_node_to_node_centroid.py " + moving_node_fname + " " + node_fname[1] + " " + "temp.node")
		os.system("mv temp.node " + moving_node_fname)
	
	
	# Run an FFEA minimisation
	os.system("ffea " + initial_script_fname)

	## Get final frame of trajectory and make a new node file for next run
	##final_frame_to_nodes(traj_fname, moving_node_fname, scale)
	
	# Change restart to 1
	if num_runs == 0:
		change_restart(initial_script_fname)

	# Ask if acceptable
	acceptance = raw_input("Please check this level of overlap in the viewer (PATH_TO_VIEWER/FFEA_viewer " + str(initial_script_fname) + "). Is this acceptable (y/n)?:")
	if acceptance == "y" or acceptance == "Y":
		not_satisfied = False
		break

	# Delete last spring if necessary
	if num_runs != 0:
		while True:
			acceptance = raw_input("Would you like to delete the last minimisation phase (y/n)?:")
			if acceptance == "y" or acceptance == "Y":
				for pair in last_set:
					node_pairs.remove(pair)
			break

	# Ask for new springs to be selected
	last_set = []
	while True:
		acceptance = raw_input("Would you like to enter a new spring (y/n)?:")
		if acceptance == "y" or acceptance == "Y":
			new_spring = [0,0]
			try:
				new_spring[0] = int(raw_input("Enter blob 0 node:"))
				new_spring[1] = int(raw_input("Enter blob 1 node:"))

			except(ValueError):
				print "Please type an integer value followed by return\n"
				continue

			node_pairs.append(new_spring)
			last_set.append(new_spring)
			num_pairs += 1
		else:
			break

	# Add simulation time
	add_sim_time(initial_script_fname, num_steps)
	make_a_spring_file(spring_fname, node_pairs, k, l)
	num_runs += 1

# Satisfied! Now make two node files and create two maps!
final_node_fname = ["final0.node", "final1.node"]
make_nodes_from_traj(final_node_fname, node_fname, traj_fname, 1000)
os.system("../Node_tools/make_structure_map " + " " + final_node_fname[0] + " " + top_fname[0] + " " + final_node_fname[1] + " " + top_fname[1] + " " + "01.map")
os.system("../Node_tools/make_structure_map " + " " + final_node_fname[1] + " " + top_fname[1] + " " + final_node_fname[0] + " " + top_fname[0] + " " + "10.map")

# Delete all new files
os.system("rm " + moving_node_fname + " " + initial_script_fname + " " + spring_fname)
