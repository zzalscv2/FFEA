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

import os, sys

# Get the path to ffeatools
FFEA_TOOLS_PATH = os.path.abspath(os.path.dirname(sys.argv[0]))
FFEA_TOOLS_PATH += "/"

# Make sure it is in the "ffeatools/" folder. If not, then it has erroneously been moved and the relative paths in the rest of this script will not work.
# if FFEA_TOOLS_PATH.endswith("/ffeatools"):
	# FFEA_TOOLS_PATH = FFEA_TOOLS_PATH + "/" # so we get the base path to ffeatools
# else:
	# sys.exit("Error: ffeatools.py must be in the ffeatools/ folder. Relative paths will not work otherwise.\n")

# Build the lookup dictionary for all FFEA tools
ffea_tools = 	{
		"pdbtoemmap": "FFEA_initialise/EM_density_map_tools/convert_pdb_to_emdb_map",
		"emmaptosurf": "FFEA_initialise/Surface_tools/Surface_convert_from_EM_density_map/emdb_map_to_ffea",
		"objtosurf": "FFEA_initialise/Surface_tools/Netgen_surf_convert_from_obj.py",
		"surftoobj": "FFEA_initialise/Surface_tools/Obj_convert_from_netgen_surf.py",
		"surftocgsurf": "FFEA_initialise/Surface_tools/surface_coarse_grainer/surface_coarse_grainer_final",
		"voltoffea": "FFEA_initialise/FFEA_convert_from_volume/FFEA_convert_from_volumetric_mesh.py",
      "cullvol": "FFEA_initialise/FFEA_volume_tools/cull_small_interior_elements.py",
		"makekineticmaps": "FFEA_initialise/FFEA_mapping_tools/FFEA_generate_kinetic_maps.py",
		"split": "FFEA_analysis/FFEA_traj_tools/FFEA_split_trajectory.py",
		"thin": "FFEA_analysis/FFEA_thin_system.py",
      "nodesFromTraj": "FFEA_analysis/FFEA_traj_tools/FFEA_get_snapshots_in_nodes.py",
      "tettonet": "FFEA_initialise/FFEA_volume_tools/convert_tet_to_net.py",
      "makestructuremap": "FFEA_initialise/FFEA_mapping_tools/make_structure_map",
      "maptosparse": "FFEA_initialise/FFEA_mapping_tools/FFEA_convert_kinetic_map_to_sparse.py",
      "maptraj": "FFEA_initialise/FFEA_mapping_tools/Kinetic_FFEA_map_apply_to_FFEA_traj.py",
      "automodel": "FFEA_initialise/FFEA_automodel.py",
      "makepseudopdb": "FFEA_analysis/FFEA_pyPca/FFEA_convert_to_PCAsystem.py",
      "pyPCAanim": "FFEA_analysis/FFEA_pyPca/FFEA_get_PCA_animations.py",
      "pyPCAeigen": "FFEA_analysis/FFEA_pyPca/FFEA_get_PCA_eigensystem.py",
      "pyPCAproj": "FFEA_analysis/FFEA_pyPca/FFEA_get_PCA_projections.py",
      "pintovdw": "FFEA_initialise/FFEA_vdw_tools/FFEA_pin_to_vdw.py"
		}

if len(sys.argv) == 1:
	usage_string = "\nUsage: ./ffeatools ACTION ARGS\n"
	usage_string += "With ACTION as one of:\n\n"
	for script in sorted(ffea_tools.keys()):
		usage_string += script+"\n"
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
