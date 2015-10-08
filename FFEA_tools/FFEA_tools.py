#!/usr/bin/env python
import os, sys

# Get the path to FFEA_tools
FFEA_TOOLS_PATH = os.path.abspath(os.path.dirname(sys.argv[0]))

# Make sure it is in the "FFEA_tools/" folder. If not, then it has erroneously been moved and the relative paths in the rest of this script will not work.
if FFEA_TOOLS_PATH.endswith("/FFEA_tools"):
	FFEA_TOOLS_PATH = FFEA_TOOLS_PATH + "/" # so we get the base path to FFEA_tools
else:
	sys.exit("Error: FFEA_tools.py must be in the FFEA_tools/ folder. Relative paths will not work otherwise.\n")

# Build the lookup dictionary for all FFEA tools
ffea_tools = 	{
		"makeffeablob": "FFEA_initialise/FFEA_convert_from_volume/FFEA_convert_from_volumetric_mesh.py",
		"makekineticmaps": "FFEA_initialise/Kinetic_FFEA_convert_from_FFEA/Kinetic_FFEA_create_maps.py",
		"meshmap": "FFEA_initialise/Surface_tools/Surface_convert_from_EM_density_map/emdb_map_to_ffea",
		"addmap": "FFEA_initialise/EM_density_map_tools/add_maps",
		"pdbtomap": "FFEA_initialise/EM_density_map_tools/convert_pdb_to_emdb_map",
		"align": "FFEA_initialise/PDB_tools/PDB_align/align_point_cloud",
		"extractnodes": "FFEA_analysis/FFEA_traj_tools/FFEA_extract_specific_nodes",
		"trajtoxyz": "FFEA_analysis/FFEA_traj_tools/FFEA_extract_positions",
		"xyztodx": "FFEA_analysis/FFEA_traj_tools/Dx_convert_from_xyz_position_trajectory.py",
		"cgsurf": "FFEA_initialise/Surface_tools/surface_coarse_grainer/surface_coarse_grainer_final",
		"cuboid": "FFEA_initialise/Volume_tools/make_cuboid_mesh/make_cuboid_mesh"
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
