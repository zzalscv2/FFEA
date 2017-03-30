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
import FFEA_vdw, FFEA_pin, FFEA_surface
import numpy as np
import __builtin__

import argparse as _argparse

# Set up argparse
parser = _argparse.ArgumentParser(description="Activate the vdw faces based on predefined pinned nodes)")
parser.add_argument("-i", action="store", nargs=1, help="Input VdW file (.vdw).")
parser.add_argument("-p", action="store", nargs=1, help="Input Pin file (.pin).")
parser.add_argument("-s", action="store", nargs=1, help="Input Surf file (.surf).")
parser.add_argument("-o", action="store", nargs=1, help="Output VdW file (.vdw).")
parser.add_argument("-ind", action="store", help="VdW Face Index")

def pin_to_vdw(vdw_fname, pin_fname, surf_fname, output_fname, index):

	# Check for problems
	if vdw_fname == None or surf_fname == None or surf_fname == None or output_fname == None:
		raise IOError

	if index == None:
		index = 0
	else:
		index = int(index)

	# Get FFEA objects
	surf = FFEA_surface.FFEA_surface(surf_fname)
	pin = FFEA_pin.FFEA_pin(pin_fname)
	vdw = FFEA_vdw.FFEA_vdw(vdw_fname)

	for i in range(surf.num_faces):
		for n in surf.face[i].n[0:3]:
			if n in pin.index:
				vdw.set_index(i, index)
				break

	vdw.write_to_file(output_fname)

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    try:
	pin_to_vdw(args.i[0], args.p[0], args.s[0], args.o[0], args.ind)
    except IOError:
	parser.print_help()
    except TypeError:
	parser.print_help()
