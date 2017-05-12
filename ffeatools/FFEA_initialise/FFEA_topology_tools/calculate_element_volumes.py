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
import FFEA_topology, FFEA_node
import numpy as np
from matplotlib import pyplot as plt
import __builtin__

import argparse as _argparse

# Set up argparse
parser = _argparse.ArgumentParser(description="Print elements in order of volume and plot")
parser.add_argument("-i", action="store", nargs=1, help="Input topology file (.top).")
parser.add_argument("-n", action="store", nargs=1, help="Input Node file (.node).")
parser.add_argument("-o", action="store", nargs=1, help="Output fname (.dat).")

def calculate_element_volumes(top_fname, node_fname, out_fname):

	# Check for problems
	if top_fname == None or node_fname == None:
		raise IOError

	# Get FFEA objects
	node = FFEA_node.FFEA_node(node_fname)
	top = FFEA_topology.FFEA_topology(top_fname)

	# Build tuples of indices and volumes
	index = 0
	tups = []
	for e in top.element:
		tups.append((index, e.calc_volume(node)))
		index += 1

	# Sort
	sorted_tups = sorted(tups, key=lambda tup: tup[1])

	# Write
	if out_fname == None:
		out_fname = os.path.splitext(top_fname)[0] + ".dat"
	if os.path.exists(out_fname):
		print("Output file already exists\n")
		raise IOError

	with open(out_fname, "w") as fout:
		fout.write("Index\tVolume\n\n")
		x = []
		y = []
		for t in sorted_tups:
			x.append(t[0])
			y.append(t[1])
			fout.write("%d\t%e\n" % (t[0], t[1]))

	# Plot
	plt.plot(y)
	plt.ylabel("Element Volumes")
	plt.title("Element Volume Range")
	plt.savefig(os.path.splitext(out_fname)[0] + ".png") 

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    try:
	calculate_element_volumes(args.i[0], args.n[0], args.o[0])
    except IOError:
	parser.print_help()
    except TypeError:
	parser.print_help()
