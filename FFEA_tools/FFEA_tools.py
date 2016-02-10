#!/usr/bin/env python
import os, sys

# Get the path to FFEA_tools
FFEA_TOOLS_PATH = os.path.abspath(os.path.dirname(sys.argv[0]))
FFEA_TOOLS_PATH += "/"

# Make sure it is in the "FFEA_tools/" folder. If not, then it has erroneously been moved and the relative paths in the rest of this script will not work.
# if FFEA_TOOLS_PATH.endswith("/FFEA_tools"):
	# FFEA_TOOLS_PATH = FFEA_TOOLS_PATH + "/" # so we get the base path to FFEA_tools
# else:
	# sys.exit("Error: FFEA_tools.py must be in the FFEA_tools/ folder. Relative paths will not work otherwise.\n")

# Build the lookup dictionary for all FFEA tools
ffea_tools = 	{
		"pdbtoemmap": "FFEA_initialise/EM_density_map_tools/convert_pdb_to_emdb_map",
		"emmaptosurf": "FFEA_initialise/Surface_tools/Surface_convert_from_EM_density_map/emdb_map_to_ffea",
		"objtosurf": "FFEA_initialise/Surface_tools/Netgen_surf_convert_from_obj.py",
		"surftoobj": "FFEA_initialise/Surface_tools/Obj_convert_from_netgen_surf.py",
		"surftocgsurf": "FFEA_initialise/Surface_tools/surface_coarse_grainer/surface_coarse_grainer_final",
		"voltoffea": "FFEA_initialise/FFEA_convert_from_volume/FFEA_convert_from_volumetric_mesh.py",
		"makekineticmaps": "FFEA_initialise/FFEA_mapping_tools/FFEA_generate_kinetic_maps.py",
		"split": "FFEA_analysis/FFEA_traj_tools/FFEA_split_trajectory.py",
		"thin": "FFEA_analysis/FFEA_traj_tools/FFEA_thin_trajectory.py",
                "view": "FFEA_analysis/FFEA_viewer/FFEA_viewer.py",
		}

if len(sys.argv) == 1:
	usage_string = "\nUsage: ./FFEA_tools ACTION ARGS\n"
	usage_string += "With ACTION as one of:\n\t"
	usage_string += " ".join(ffea_tools.keys()) + "\n"
	sys.exit(usage_string)

action = sys.argv[1]
args = " ".join(sys.argv[2:])
if action not in ffea_tools.keys():
	sys.exit("Unrecognised ACTION '" + action + "'\nACTION must be one of:\n\t" + " ".join(ffea_tools.keys()) + "\n")

# if it is a python script: 
if ffea_tools[action][-3:] == ".py":
  command = "python " + FFEA_TOOLS_PATH + ffea_tools[action] + " " + args
else:
  command = FFEA_TOOLS_PATH + ffea_tools[action] + " " + args

print "ACTION = '" + action + "' => Using command '" + command + "'\n"

os.system(command + "\n")
