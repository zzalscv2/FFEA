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
parser = _argparse.ArgumentParser(description="Activate the vdw faces based on distance from a predefined face)")
parser.add_argument("-i", action="store", nargs=1, help="Input VdW file (.vdw).")
parser.add_argument("-n", action="store", nargs=1, help="Input Node file (.node).")
parser.add_argument("-s", action="store", nargs=1, help="Input Surf file (.surf).")
parser.add_argument("-o", action="store", nargs=1, help="Output VdW file (.vdw).")
parser.add_argument("-r", action="store", help="Radius of activation")
parser.add_argument("-find", action="store", help="Central Face Index")
parser.add_argument("-ind", action="store", help="VdW Face Type (-1 -> 7)")

def activation_radially_from_face(vdw_fname, node_fname, surf_fname, output_fname, rad, findex, index):

	# Check for problems
	if vdw_fname == None or node_fname == None or surf_fname == None or output_fname == None or rad == None or findex == None:
		raise IOError

	if index == None:
		index = 0
	else:
		index = int(index)

	findex = int(findex)

	# Get FFEA objects
	surf = FFEA_surface.FFEA_surface(surf_fname)
	node = FFEA_node.FFEA_node(node_fname)
	vdw = FFEA_vdw.FFEA_vdw(vdw_fname)

	# Get distance of centroid to core node
	rad = float(rad)
	fcent = surf.face[findex].calc_centroid(node)
	for i in range(surf.num_faces):
		c = surf.face[i].calc_centroid(node)
		if np.linalg.norm(c - fcent) <= rad:
			vdw.set_index(i, index)

	vdw.write_to_file(output_fname)

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    try:
        activation_radially_from_face(args.i[0], args.n[0], args.s[0], args.o[0], args.r, args.find, args.ind)
    except IOError:
        parser.print_help()
    except TypeError:
        parser.print_help()
	raise IOError("Likely missing argument. Please try again :)")
