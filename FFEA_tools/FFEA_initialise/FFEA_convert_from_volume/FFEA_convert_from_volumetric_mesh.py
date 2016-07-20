#!/usr/bin/env python

import os, sys
from FFEA_universe import *

def get_num_nodes(nodes_fname):
	with open(nodes_fname, "r") as f:
		for line in f:
			line = line.split()
			if "num_nodes" in line[0]:
				return int(line[1])

def get_num_elements(top_fname):
	with open(top_fname, "r") as f:
		for line in f:
			line = line.split()
			if "num_elements" in line[0]:
				return int(line[1])

inputfile = None
node = None
top = None
surf = None
scaling_factor = "1.0"
stokes = None
stokes_radius = None
material = None
mat_rho = None
mat_shear_visc = None
mat_bulk_visc = None
mat_shear_mod = None
mat_bulk_mod = None
mat_dielec = None
quality = None
small_cull = None
make_script = False

no_stokes = True
calc_stokes = False
no_mat = True
i = 1

while i < len(sys.argv):
	if sys.argv[i] == "-mesh":
		i += 1
		inputfile = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-node":
		i += 1
		node = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-top":
		i += 1
		top = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-surf":
		i += 1
		surf = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-scale":
		i += 1
		scaling_factor = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-stokes":
		i += 1
		stokes = sys.argv[i]
		no_stokes = False
		i += 1
		continue
	if sys.argv[i] == "-stokes_radius":
		i += 1
		stokes_radius = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-material":
		i += 1
		material = sys.argv[i]
		no_mat = False
		i += 1
		continue
	if sys.argv[i] == "-density":
		i += 1
		mat_rho = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-shear_visc":
		i += 1
		mat_shear_visc = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-bulk_visc":
		i += 1
		mat_bulk_visc = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-shear_mod":
		i += 1
		mat_shear_mod = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-bulk_mod":
		i += 1
		mat_bulk_mod = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-dielec":
		i += 1
		mat_dielec = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-quality":
		i += 1
		quality = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-small_cull":
		i += 1
		small_cull = sys.argv[i]
		i += 1
		continue
	if sys.argv[i] == "-make_script":
		make_script = True
		i += 1
		continue
	sys.exit("Unrecognised option '"+ sys.argv[i] +"'")

if inputfile == None:
	sys.exit("Please specify an inputfile (.vol, .off or .dat) using -mesh")

if stokes == None and stokes_radius == None:
	if(raw_input("No stokes file or homogeneous stokes radius specified.\nWould you like us to calculate an optimum radius (y or n)?") == "y"):
		calc_stokes = True
		no_stokes = False
	else:
		sys.exit("Please specify a stokes file using -stokes, or a homogeneous stokes radius using -stokes_radius")

if material == None and (mat_rho == None or mat_shear_visc == None or mat_bulk_visc == None or mat_shear_mod == None or mat_bulk_mod == None or mat_dielec == None):
	sys.exit("Please specify a material file using -material, or specify ALL homogeneous parameters using -density, -shear_visc, -bulk_visc, -shear_mod, -bulk_mod and -dielec")

path = os.path.dirname(sys.argv[0])
if path == '':
	path = '.'

print "path of running script: " + path

print "Input file = " + inputfile
root, ext = os.path.splitext(inputfile)
base = os.path.basename(root)
print "root = " + root
print "extension = " + ext
print "base = " + base

if node == None:
	node = base + ".node"	
if top == None:
	top = base + ".top"	
if surf == None:
	surf = base + ".surf"	
if stokes == None:
	stokes = base + ".stokes"	
if material == None:
	material = base + ".mat"	

vdw = base + ".vdw"
ljmat = base + ".lj"
default_script_fname = base + ".ffea"
pin = base + ".pin"

print "Script will run with the following parameters:"
print "node file = " + node
print "topology file = " + top
print "surface file = " + surf
print "stokes file = " + stokes
print "material file = " + material
print "scaling factor = " + scaling_factor

print "Determining file type from extension:"
if ext == ".vol":

	# Check consistency of file
	topology = FFEA_topology.FFEA_topology(inputfile)
	if topology == None or topology.num_elements == 0:
		sys.exit("Error. No elements found within '" + inputfile + "'. Please provide a consistent file.")

	if small_cull == None:
		print "VOL extension = Netgen mesh output file. Converting to walrus format..."
		os.system("python " + path + "/netgen_to_walrus.py " + inputfile + " " + node + " " + top + "\n")
	else:
		os.system("python " + path + "/convert_vol_to_tetgen_output.py " + inputfile + " " + root + ".1.node " + root + ".1.ele\n")
		print "Culling all elements below volume " + small_cull
		os.system("python " + path + "/cull_small_elements.py " + root + ".1.node " + root + ".1.ele " + small_cull + " " + scaling_factor + "\n")
		print "Now converting tetgen output to walrus format..."
		os.system("python " + path + "/tetgen_to_walrus.py " + root + ".1.node " + root + ".1.ele " + node + " " + top + "\n")
elif ext == ".off":
	print "OFF extension = DEC Object File Format."
	print "Warning: This is a surface mesh, and needs to be meshed."
	if quality == None:
		print "Meshing with 'tetgen -a -Q " + inputfile + "'..."
		os.system("tetgen -a -Q " + inputfile)
	else:
		print "Meshing with 'tetgen -e " + inputfile + "'..."
		os.system("tetgen -e " + inputfile)
	if small_cull != None:
		print "Culling all elements below volume " + small_cull
		os.system("python " + path + "/cull_small_elements.py " + root + ".1.node " + root + ".1.ele " + small_cull + " " + scaling_factor + "\n")
	print "Now converting tetgen output to walrus format..."
	os.system("python " + path + "/tetgen_to_walrus.py " + root + ".1.node " + root + ".1.ele " + node + " " + top + "\n")
elif ext == ".dat":
	print "DAT extension = Mesh format"
	os.system("python " + path + "/DAT_to_walrus.py " + inputfile + " " + node + " " + top + "\n")
else:
	print "Error: File extension not recognised. Only '.vol' and '.off' are supported."
	sys.exit(0)


# create stokes file if necessary
if calc_stokes == True:
	print "Calculating homogeneous stokes file..."
	os.system("python " + path + "/calc_homo_stokes_file.py " + node + " " + stokes + "\n")
if no_stokes == True:
	print "Creating homogeneous stokes file..."
	num_nodes = get_num_nodes(node)
	os.system("python " + path + "/create_homo_stokes_file.py " + str(num_nodes) + " " + stokes_radius + " " + stokes + "\n")

# create material file if necessary
if no_mat == True:
	num_elements = get_num_elements(top)
	os.system("python " + path + "/create_homo_material_params_file.py " + str(num_elements) + " " + mat_rho + " " + mat_shear_visc + " " + mat_bulk_visc + " " + mat_shear_mod + " " + mat_bulk_mod + " " + mat_dielec + " " + material + "\n")


# convert 4 node tets to 10 point tets
os.system("python " + path + "/convert_tetrahedra_linear_to_quadratic.py " + node + " " + top + " " + stokes + "\n")
#sys.exit()

# extract the surface
os.system("python " + path + "/extract_surface_from_topology.py " + top + " " + surf + "\n")
#sys.exit()
# move all surface nodes to the top of the list
os.system("python " + path + "/put_surface_nodes_first.py " + node + " " + top + " " + surf + " " + stokes + "\n")
#sys.exit()
# move all surface elements to the top of the list
os.system("python " + path + "/put_surface_elements_first.py " + node + " " + top + " " + surf + " " + material + "\n")

# cycle surface indices so that all normals point out of the mesh
os.system("python " + path + "/make_all_normals_point_outwards.py " + node + " " + top + " " + surf + "\n")

# calculate the element connectivity (for use with van der waals calculations)
#os.system("python " + path + "/calculate_element_connectivity.py " + top + " " + vdw + "\n")

# scale the cube by the desired factor
os.system("python " + path + "/scale.py " + node + " " + node + " " + scaling_factor + "\n")

# make a default noninteracting vdw file
os.system("python " + path + "/make_default_vdw.py " + surf + " " + vdw + "\n")

# make a default vdw interaction matrix file
os.system("python " + path + "/make_default_LJmatrix.py " + ljmat + "\n")

# make a default no pinning pin file
os.system("python " + path + "/make_default_pin.py " + pin + "\n")

# if requested, make a default script for this one blob
if make_script == True:
	os.system("python " + path + "/make_default_script.py " + base + " " + default_script_fname + "\n")

print "All done"
