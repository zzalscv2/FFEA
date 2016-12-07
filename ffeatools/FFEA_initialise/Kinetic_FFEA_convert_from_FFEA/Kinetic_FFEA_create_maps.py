#!/usr/bin/env python
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


import sys, os

# Function for getting rid of all files created by this program
def delete_files():
	
	os.system("rm " + working_spring_fname + " " + working_script_fname + " " + "working.pin" + " " + "final0.node" + " " + "final1.node")

# Function for getting files from script
def get_structure_from_script(fname):
	
	# Open script and get lines
	fin = open(fname, "r")
	lines = fin.readlines()
	fin.close()

	# Get structure files
	node_fname = []
	top_fname = []
	for line in lines:
		if "=" not in line:
			continue
		lvalue, rvalue = line.strip()[1:-1].split("=")
		if lvalue.strip() == "nodes":
			node_fname.append(rvalue.strip())
		elif lvalue.strip() == "topology":
			top_fname.append(rvalue.strip())

	# Quit if not enough files
	if len(node_fname) != 2 or len(top_fname) != 2:
		delete_files()
		sys.exit("Error. '" + fname + "' does not contain two blobs as required.")

	return node_fname[0], node_fname[1], top_fname[0], top_fname[1]

def get_parameters_from_script(fname, param_lvalue):
	
	# Open script and get lines
	fin = open(fname, "r")
	lines = fin.readlines()
	fin.close()
	
	# Get values
	parameters = []
	for line in lines:
		if "=" not in line:
			continue
		lvalue, rvalue = line.strip()[1:-1].split("=")
		if lvalue.strip() == param_lvalue:
			parameters.append(rvalue.strip())
	
	return parameters

# Main program begins here
if len(sys.argv) < 5:
	sys.exit("Usage python " + sys.argv[0] + " [INPUT .ffea script (must have 1 blob w/ 2 conformations)] [OUTPUT .map fname] [spring_constant] [tolerance length] [output check rate (num_frames)] [Pairs of nodes {[i,j]}]\n")

# Get initial scripts
py_script_dirname = os.path.dirname(sys.argv[0])
script_fname = sys.argv[1]
in_script_base, in_script_ext = os.path.splitext(os.path.abspath(script_fname))

map_fname = sys.argv[2]
out_map_base, out_map_ext = os.path.splitext(os.path.abspath(map_fname))
if out_map_ext == "":
	out_map_ext = ".map"

out_map_fname = []
out_map_fname.append(out_map_base + "_01_dense" + out_map_ext)
out_map_fname.append(out_map_base + "_10_dense" + out_map_ext)
out_map_fname.append(out_map_base + "_01_sparse" + out_map_ext)
out_map_fname.append(out_map_base + "_10_sparse" + out_map_ext)

k = float(sys.argv[3])
l = float(sys.argv[4])
check = int(sys.argv[5])
node_pairs = []
node_pair_string = ""

if len(sys.argv) > 5:
	for i in range(6, len(sys.argv), 1):
		try:
			node_pairs.append([int(sys.argv[i].strip()[1:-1].split(",")[0]), int(sys.argv[i].strip()[1:-1].split(",")[1])])
			node_pair_string += sys.argv[i].strip()[1:-1].split(",")[0] + " " + sys.argv[i].strip()[1:-1].split(",")[1] + " "
		except IndexError:
			sys.exit("Error. '" + sys.argv[i] + "' improperly formatted. Node pairs should look like this: [i,j]")
		except ValueError:
			sys.exit("Error. '" + sys.argv[i] + "' contains a non integer character. Or improperly formatted. Node pairs should look like this: [i,j]")

# Get node and topology fnames
orig_node_fname = ["",""]
orig_top_fname = ["",""]
orig_node_fname[0], orig_node_fname[1], orig_top_fname[0], orig_top_fname[1] = get_structure_from_script(script_fname)

# Make spring file
if len(node_pairs) != 0:
	working_spring_fname = in_script_base + ".spring"
	while os.path.isfile(working_spring_fname):
		working_spring_fname = working_spring_fname[0:-7] + "_a" + ".spring"

	os.system("python " + py_script_dirname + "/Kinetic_FFEA_create_working_spring_file.py " + working_spring_fname + " " + str(k) + " " + str(l) + " " + node_pair_string)
else:
	working_spring_fname = "none.spring"

# Make working script
working_script_fname = in_script_base + "_working.ffea"
while os.path.isfile(working_script_fname):
	working_script_fname = working_script_fname[0:-5] + "_a" + ".ffea"
os.system("python " + py_script_dirname + "/Kinetic_FFEA_create_working_script.py " + script_fname + " " + working_script_fname + " " + working_spring_fname + " " + "0" + " " + str(check))

# Run working script repeatadly until user is satisfied
run = 0
while True:

	# Make restart after first go
	if run > 0:
		os.system("python " + py_script_dirname + "/Kinetic_FFEA_create_working_script.py " + script_fname + " " + working_script_fname + " " + working_spring_fname + " " + str(run) + " " + str(check))

	os.system("ffea " + working_script_fname)
	run += 1
	finished = raw_input("Are you satisfied with this level of overlap (y/n)? (check '" + os.path.basename(working_script_fname) + "' in the viewer):")
	if finished.lower() == "y":
		finished = raw_input("Shall we create the maps (y/n)?:")
		if finished.lower() == "y":
			break
		else:
			delete_files()
			sys.exit("Exiting without creating maps.\n")

	extra_spring = raw_input("Would you like to add another spring (y/n)?:")
	if extra_spring.lower() == "y":
		working_spring_fname = in_script_base + ".spring"
		os.system("rm none")
		while True:

			line = raw_input("Enter spring node pair, or enter to finish:")
			if line == "":
				break
			
			sline = line[1:-1].split(",")
			extra_node_pair = [int(sline[0]), int(sline[1])]
			node_pairs.append(extra_node_pair)
			node_pair_string += str(extra_node_pair[0]) + " " + str(extra_node_pair[1]) + " "

		os.system("python " + py_script_dirname + "/Kinetic_FFEA_create_working_spring_file.py " + working_spring_fname + " " + str(k) + " " + str(l) + " " + node_pair_string)
		continue

	delete_spring = raw_input("Would you like to delete the last spring (y/n)?:")
	if delete_spring.lower() == "y":
		node_pairs.pop()
		node_pair_string = node_pair_string[0:-4]
		os.system("python /localhome/py09bh/Software/FFEA/FFEA_git/ffeatools/FFEA_initialise/Kinetic_FFEA_convert_from_FFEA/Kinetic_FFEA_create_working_spring_file.py " + working_spring_fname + " " + str(k) + " " + str(l) + " " + node_pair_string)
		continue

# Satisfied! Now make two node files and create two maps!
final_node_fname = ["final0.node", "final1.node"]
working_traj_fname = os.path.abspath(get_parameters_from_script(working_script_fname, "trajectory_out_fname")[0])
scale = [1 for i in range(2)]
scale = get_parameters_from_script(working_script_fname, "scale")
os.system("python " + py_script_dirname + "/../../FFEA_analysis/FFEA_traj_tools/FFEA_traj_to_nodes.py " + working_traj_fname + " " + "1000" + " " + " 0 0 " + orig_node_fname[0] + " " + final_node_fname[0] + " " + scale[0] + " " + " 1 0 " + orig_node_fname[1] + " " + final_node_fname[1] + " " + scale[0])
os.system(py_script_dirname + "/../FFEA_node_tools/make_structure_map " + " " + final_node_fname[0] + " " + orig_top_fname[0] + " " + final_node_fname[1] + " " + orig_top_fname[1] + " " + out_map_fname[0])
os.system(py_script_dirname + "/../FFEA_node_tools/make_structure_map " + " " + final_node_fname[1] + " " + orig_top_fname[1] + " " + final_node_fname[0] + " " + orig_top_fname[0] + " " + out_map_fname[1])

# Convert maps to sparse format
os.system("python " + py_script_dirname + "/../FFEA_kinetic_map_tools/Kinetic_FFEA_map_convert_to_sparse.py " + out_map_fname[0] + " " + out_map_fname[2])
os.system("python " + py_script_dirname + "/../FFEA_kinetic_map_tools/Kinetic_FFEA_map_convert_to_sparse.py " + out_map_fname[1] + " " + out_map_fname[3])

# Remove stuff that was made with this script
delete_files()
