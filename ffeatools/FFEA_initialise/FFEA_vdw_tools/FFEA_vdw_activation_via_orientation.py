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

import sys, os
import FFEA_vdw, FFEA_node, FFEA_surface
import numpy as np
import __builtin__

import argparse as _argparse

# Set up argparse
parser = _argparse.ArgumentParser(description="Activate the vdw faces based orientation of faces)")
parser.add_argument("-i", action="store", nargs=1, help="Input VdW file (.vdw).")
parser.add_argument("-n", action="store", nargs=1, help="Input Node file (.node).")
parser.add_argument("-s", action="store", nargs=1, help="Input Surf file (.surf).")
parser.add_argument("-o", action="store", nargs=1, help="Output VdW file (.vdw).")
parser.add_argument("-vn", action="store", nargs=3, help="Orientation vector")
parser.add_argument("-tol", action="store", nargs=1, help="Tolerance (degrees)")
parser.add_argument("-ind", action="store", help="VdW Face Index")

def activate_via_orientation(vdw_fname, node_fname, surf_fname, output_fname, vn, tol, index):

	# Check for problems
	if vdw_fname == None or node_fname == None or surf_fname == None or output_fname == None:
		raise IOError
	if vn == None or len(vn) != 3:
		raise IOError

	if index == None:
		index = 0
	else:
		index = int(index)

	# Build objects
	vdw = FFEA_vdw.FFEA_vdw(vdw_fname)
	node = FFEA_node.FFEA_node(node_fname)
	surf = FFEA_surface.FFEA_surface(surf_fname)	

	try:
		# Normalise vn
		vn = np.array([float(i) for i in vn])
		vn *= 1.0 / np.linalg.norm(vn)
		tol = float(tol)

	except(ValueError):
		raise

	# For all faces, if orientation vector is within lim degrees of face normal, activate it
	for i in range(surf.num_faces):

		# Calculate side of plane
		n = surf.face[i].calc_normal(node)	
		ang = np.arccos(np.dot(n, vn)) * (180.0 / np.pi)
		if ang <= tol:
			vdw.set_index(i, index)
		
	# Output
	vdw.write_to_file(output_fname)

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    try:
	activate_via_orientation(args.i[0], args.n[0], args.s[0], args.o[0], args.vn, args.tol[0], args.ind)
    except IOError:
	parser.print_help()
    except TypeError:
	parser.print_help()
    except ValueError:
	parser.print_help()
