#!/usr/bin/env python
import os, sys

# Get the path to FFEA_tools
FFEA_TOOLS_PATH = os.path.abspath(os.path.dirname(sys.argv[0]))

# Make sure it is in the "all/" folder. If not, then it has erroneously been moved and the relative paths in the rest of this script will not work.
if FFEA_TOOLS_PATH.endswith("/all"):
	FFEA_TOOLS_PATH = FFEA_TOOLS_PATH[:-3] # chop off the end so we get the base path to FFEA_tools
else:
	sys.exit("Error: FFEA_tools.py must be in the all/ folder in FFEA_tools. Relative paths will not work otherwise.\n")

print FFEA_TOOLS_PATH

# Build the lookup dictionary for all FFEA tools
ffea_tools = 	{
		"makeffea": "convert_to_FFEA/convert_all_to_FFEA.py",
		"meshmap": "emdb_to_surface_mesh/emdb_map_to_ffea",
		"addmap": "emdb_to_surface_mesh/add_maps",
		"pdbtomap": "emdb_to_surface_mesh/convert_pdb_to_emdb_map",
		"align": "align_point_clouds/align_point_cloud",
		"extractnodes": "extract/extract",
		"trajtoxyz": "extract/ffea_traj_to_xyz",
		"xyztodx": "extract/xyz_to_dx_map.py",
		"cgsurf": "surface_coarse_grainer/surface_coarse_grainer_final",
		"comparemodes": "compare_modes/compare_modes.py",
		"cuboid": "make_cuboid_mesh/make_cuboid_mesh"
		}

if len(sys.argv) == 1:
	usage_string = "\nUsage: ./FFEA_tools ACTION ARGS\n"
	usage_string += "With ACTION as one of:\n"
	usage_string += " ".join(ffea_tools.keys()) + "\n"
	sys.exit(usage_string)

action = sys.argv[1]
args = " ".join(sys.argv[2:])
if action not in ffea_tools.keys():
	sys.exit("Unrecognised ACTION '" + action + "'\nACTION must be one of:\n" + " ".join(ffea_tools.keys()) + "\n")

command = FFEA_TOOLS_PATH + ffea_tools[action] + " " + args
print "ACTION = '" + action + "' => Using command '" + command + "'\n"

os.system(command + "\n")
