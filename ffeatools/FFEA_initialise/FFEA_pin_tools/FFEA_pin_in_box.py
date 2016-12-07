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

import sys
import FFEA_node, FFEA_pin
from numpy import linalg

if len(sys.argv) != 9:
	sys.exit("Usage: python FFEA_pin_in_box.py [INPUT .node fname] [OUTPUT .pin fname] [x limits (max min)] [y limits (max min)] [z limits (max min)]")

# Get args
node_fname = sys.argv[1]
pin_fname = sys.argv[2]
limits = [float(sys.argv[i + 3]) for i in range(6)]

# Get nodes and create pinned node object
nodes = FFEA_node.FFEA_node(node_fname)
pinned = FFEA_pin.FFEA_pin(pin_fname)
pinned.reset()

# Set nodal centroid
centroid = [0.0, 0.0, 0.0]
nodes.set_centroid(centroid)

# Pin some nodes!
for i in range(nodes.num_nodes):
	if nodes.pos[i][0] >= limits[0] and nodes.pos[i][0] <= limits[1] and nodes.pos[i][1] >= limits[2] and nodes.pos[i][1] <= limits[3] and nodes.pos[i][2] >= limits[4] and nodes.pos[i][2] <= limits[5]:
		pinned.add_node(i)

# Output
pinned.write_to_file(pin_fname)
