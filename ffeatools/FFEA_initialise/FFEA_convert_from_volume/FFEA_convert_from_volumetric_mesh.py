#
# This script converts the .vol file from the initialisation routines into the necessary file formats for an FFEA simulation
# This means get linear element from vol and make 2nd order, store new faces, tets etc, build vdw, pin, bsites, stokes, move blob to centroid
#

# Node - 2nd order
# Top - 2nd order
# Surf - 2nd order
# Mat - 1st order
# Pin - N/A or 2nd order (linked to nodes)
# VdW - 2nd Order
# lj - N/A
# Stokes - 2nd order
# Bsites - N/A or 2nd order (linked to surf)
# beads - N/A

def error():
	sys.stdout.write("\033[91mERROR:\033[0m ")

def print_usage(av):
	print "\nUsage: python " + os.path.basename(os.path.abspath(av[0])) + " --mesh [INPUT .vol] --out [OUTPUT .ffea]"

def print_help(av):
	print_usage(av)
	print "\nOptions:\n"
	help = {"help": "Prints this help message", "mesh": "input .vol file", "stokes": "a stokes radius for each node", "cull": "removes interior elements with volume less than that specified"}
	matparams = {"density": "Density", "shear-viscosity": "Shear Viscosity", "bulk-viscosity": "Bulk Viscosity", "shear-modulus": "Shear Modulus", "bulk-modulus": "Bulk Modulus", "dielectric": "Dielectric Constant"}

	for key in help:
		print "\t--" + str(key) + "\t\t\t" + help[key]

	print "\nMaterial Parameters\n"
	for key in matparams:
		print "\t--" + str(key) + "\t\t" + matparams[key]

	sys.exit("\n")


import sys, os
import FFEA_script
from FFEA_universe import *

av = sys.argv

if len(av) < 2:
	print_usage(av)
	sys.exit()

# Default args
volfname = ""
outfname = ""
make_script = False
matparams = {"d": 1.5e3, "sv": 1e-3, "bv": 1e-3, "sm": 370e6, "bm": 111e7, "di": 1.0}
stokes_radius = None
cull = [False, 0.0]

# Get args
i = 1
while i < len(av):
	if av[i] == "--mesh" or av[i] == "-m":
		volfname = av[i + 1]
	elif av[i] == "--out" or av[i] == "-o":
		outfname = av[i + 1]
	elif av[i] == "--make-script":
		make_script = True
		i -= 1

	elif av[i] == "--help" or av[i] == "-h":
		print_help(av)
	elif av[i] == "--density" or av[i] == "-d":
		matparams["d"] = float(av[i + 1])
	elif av[i] == "--shear-viscosity" or av[i] == "--shear-visc" or av[i] == "-sv":
		matparams["sv"] = float(av[i + 1])
	elif av[i] == "--bulk-viscosity" or av[i] == "--bulk-visc" or av[i] == "-bv":
		matparams["bv"] = float(av[i + 1])
	elif av[i] == "--shear-modulus" or av[i] == "--shear-mod" or av[i] == "-sm":
		matparams["sm"] = float(av[i + 1])
	elif av[i] == "--bulk-modulus" or av[i] == "--bulk-mod" or av[i] == "-bm":
		matparams["bm"] = float(av[i + 1])
	elif av[i] == "--dielectric" or av[i] == "--dielec" or av[i] == "-di":
		matparams["di"] = float(av[i + 1])
	elif av[i] == "--stokes" or av[i] == "--stokes-radius" or av[i] == "-s":
		stokes_radius = float(av[i + 1])
	elif av[i] == "--cull" or av[i] == "-c":
		cull[0] = True
		cull[1] = float(av[i + 1])
	else:
		pass

	i += 2

# Test args
if volfname == "":
	error()
	sys.exit("You must specify a .vol filename with '--mesh'")

if outfname == "":
	basename = os.path.splitext(os.path.abspath(volfname))[0]
else:
	basename = os.path.splitext(os.path.abspath(outfname))[0]
	
# Get initial stuff from .vol file!
node = FFEA_node.FFEA_node(volfname)
top = FFEA_topology.FFEA_topology(volfname)
surf = FFEA_surface.FFEA_surface(volfname)

# Let each surface face know which element it is connected to (if this is slow, just load surface from topology instead)
if surf.get_element_indices(top) == -1:
	surf = top.extract_surface()

# Now, build necessary things that only have linear properties
mat = FFEA_material.FFEA_material()
mat.build(top.num_elements, d=matparams["d"], sv=matparams["sv"], bv=matparams["bv"], sm=matparams["sm"], bm=matparams["bm"], di=matparams["di"])

# Convert stuff to 2nd order
top.increase_order(node=node, surf=surf)

# Find out what is interior and what is surface and reorder stuff
top.calculateInterior(surf=surf)
node.calculateInterior(top=top, surf=surf)

# Check the normals in the and surface files only (this is all that is necessary for FFEA)
surf.check_normals(node, top)

# Cull small elements
if cull[0]:
	top.cull_interior(cull[1], node, surf=surf)

# Everything should be set! Make default files for the other stuff
pin = FFEA_pin.FFEA_pin()
vdw = FFEA_vdw.FFEA_vdw()
vdw.default(surf.num_faces)
lj = FFEA_lj.FFEA_lj()
lj.default()

# User may set a stokes radius. If they don't, make a default one
if stokes_radius == None:

	# All drag is currently local, so get largest length scale and treat as sphere. Very bad approximation for long molecules
	dims = node.calculate_dimensions()
	stokes_radius = max(dims) / len(top.get_linear_nodes())

stokes = FFEA_stokes.FFEA_stokes()
stokes.default(node.num_nodes, top, stokes_radius)

# Script will be defaulted
if make_script:
	script = FFEA_script.FFEA_script()
	script.default(basename)

# Now, print them all out!
node.write_to_file(basename + ".node")
top.write_to_file(basename + ".top")
surf.write_to_file(basename + ".surf")
mat.write_to_file(basename + ".mat")
vdw.write_to_file(basename + ".vdw")
lj.write_to_file(basename + ".lj")
pin.write_to_file(basename + ".pin")
stokes.write_to_file(basename + ".stokes")

if make_script:
	script.write_to_file(basename + ".ffea")
