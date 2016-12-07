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

if len(sys.argv) < 2:
	sys.exit("Usage: python " + sys.argv[0] + " [OUTPUT .pin fname] [List of nodes to pin]")

outpinfname = sys.argv[1]
nodes_to_pin = []

for i in range(2, len(sys.argv)):
	nodes_to_pin.append(sys.argv[i])

# Pin all nodes
fout = open(outpinfname, "w")
fout.write("ffea pinned nodes file\n")
fout.write("num_pinned_nodes %d\n" % (len(nodes_to_pin)))
fout.write("pinned nodes:\n")
for i in range(len(nodes_to_pin)):
	fout.write(nodes_to_pin[i] + "\n")
fout.close()	
