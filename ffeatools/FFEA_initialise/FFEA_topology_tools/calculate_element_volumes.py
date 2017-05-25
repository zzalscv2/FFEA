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
		tups.append((index, e.calc_volume(node), e.get_smallest_lengthscale(node)))
		index += 1

	# Sort
	sorted_tupsvol = sorted(tups, key=lambda tup: tup[1])
	sorted_tupslen = sorted(tups, key=lambda tup: tup[2])

	# Write
	if out_fname == None:
		out_fname = os.path.splitext(top_fname)[0] + ".dat"
#	if os.path.exists(out_fname):
#		print("Output file already exists\n")
#		raise IOError

	with open(out_fname, "w") as fout:
		fout.write("Index\tVolume\t\t\tIndex\tLength\n\n")
		x1 = []
		y1 = []
		x2 = []
		y2 = []
		for t1,t2 in zip(sorted_tupsvol, sorted_tupslen):
			x1.append(t1[0])
			y1.append(t1[1])
			x2.append(t2[0])
			y2.append(t2[2])
			fout.write("%d\t%e\t\t%d\t%e\n" % (t1[0], t1[1], t2[0], t2[2]))

	# Write important stuff
	print("\nElement Volume Details for '" + top_fname + "':\n")
	print("\tSmallest: Index=%d, Volume=%f, Length=%f" % (x1[0], y1[0], top.element[x1[0]].get_smallest_lengthscale(node)))
	print("\tLargest: Index=%d, Volume=%f, Length=%f" % (x1[-1], y1[-1], top.element[x1[-1]].get_smallest_lengthscale(node)))
	print("\tAverage: %f +/- %f" % (np.mean(y1), np.std(y1)))
	
	print("\n\nElement Length Details for '" + top_fname + "':\n")
	print("\tSmallest: Index=%d, Length=%f, Volume=%f" % (x2[0], y2[0], top.element[x2[0]].calc_volume(node)))
	print("\tLargest: Index=%d, Volume=%f, Length=%f" % (x2[-1], y2[-1], top.element[x2[-1]].calc_volume(node)))
	print("\tAverage: %f +/- %f" % (np.mean(y2), np.std(y2)))

	# Plot
	plt.figure()
	plt.plot(y1)
	plt.ylabel("Element Volumes")
	plt.title("Element Volume Range")
	plt.savefig(os.path.splitext(out_fname)[0] + "_vols.png")

	plt.figure()
	plt.plot(y2)
	plt.ylabel("Element Lengths")
	plt.title("Element Length Range")
	plt.savefig(os.path.splitext(out_fname)[0] + "_lengths.png") 

if sys.stdin.isatty() and hasattr(__builtin__, 'FFEA_API_mode') == False:
    args = parser.parse_args()
    try:
        calculate_element_volumes(args.i[0], args.n[0], args.o[0])
    except IOError:
	parser.print_help()
    except TypeError:
	parser.print_help()
