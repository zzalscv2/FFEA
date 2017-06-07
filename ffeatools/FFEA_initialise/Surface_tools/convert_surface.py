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
import FFEA_surface, FFEA_node
import __builtin__

import argparse as _argparse

# Set up argparse
parser = _argparse.ArgumentParser(description="Convert between surface formats depending on file extensions)")
parser.add_argument("-i", action="store", nargs=1, help="Input Surface file.")
parser.add_argument("-o", action="store", nargs=1, help="Input Surface file.")

def convert_surface(insfname, outsfname):

	# Check for problems
	if insfname == None or outsfname == None:
		raise IOError
	
	# Try to build and save based on the file indices. Objects should raise appropriate errors if something is wrong
	surf = FFEA_surface.FFEA_surface(insfname)
	node = FFEA_node.FFEA_node(insfname)
	surf.write_to_file(outsfname, node=node)

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    try:
	convert_surface(args.i[0], args.o[0])
    except IOError:
	parser.print_help()
    except TypeError:
	parser.print_help()
