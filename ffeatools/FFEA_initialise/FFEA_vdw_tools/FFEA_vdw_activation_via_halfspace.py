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
import FFEA_vdw, FFEA_node, FFEA_surface
import numpy as np
import __builtin__

import argparse as _argparse

# Set up argparse
parser = _argparse.ArgumentParser(description="Activate the vdw faces based on one side of a defined plane)")
parser.add_argument("-i", action="store", nargs=1, help="Input VdW file (.vdw).")
parser.add_argument("-n", action="store", nargs=1, help="Input Node file (.node).")
parser.add_argument("-s", action="store", nargs=1, help="Input Surf file (.surf).")
parser.add_argument("-o", action="store", nargs=1, help="Output Node file (.vdw).")
parser.add_argument("-vn", action="store", nargs=3, help="Plane normal vector")
parser.add_argument("-vp", action="store", nargs=3, help="Position in plane")
parser.add_argument("-ind", action="store", help="VdW Face Index")

def activate_via_halfspace(vdw_fname, node_fname, surf_fname, output_fname, vn, vp, index):

	# Check for problems
	if vdw_fname == None or node_fname == None or surf_fname == None or output_fname == None:
		raise IOError
	if vn == None or len(vn) != 3:
		raise IOError
	if vp == None or len(vp) != 3:
		raise IOError

	if index == None:
		index = 0
	else:
		index = int(index)

	# Build objects
	vdw = FFEA_vdw.FFEA_vdw(vdw_fname)
	node = FFEA_node.FFEA_node(node_fname)
	surf = FFEA_surface.FFEA_surface(surf_fname)	

	vn = np.array([float(i) for i in vn])
	vp = np.array([float(i) for i in vp])

	# For each face
	for i in range(surf.num_faces):
		
		# Calculate side of plane
		c = surf.face[i].calc_centroid(node)
		
		pton = c - vp

		if np.dot(pton, vn) >= 0:
			vdw.set_index(i, index)

	# Output
	vdw.write_to_file(output_fname)

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    try:
	activate_via_halfspace(args.i[0], args.n[0], args.s[0], args.o[0], args.vn, args.vp, args.ind)
    except IOError:
	parser.print_help()
    #convert_from_volumetric_mesh(args.mesh, args.stokes_radius, args.cull, args.density, args.shear_visc, args.bulk_visc, args.shear_mod, args.bulk_mod, args.dielec, args.make_script, args.out)
