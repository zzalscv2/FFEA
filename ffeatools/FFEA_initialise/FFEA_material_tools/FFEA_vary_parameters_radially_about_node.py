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

import FFEA_material, FFEA_node, FFEA_topology
import sys, os
import numpy as np

if len(sys.argv) != 13:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT .mat file] [OUTPUT .mat file] [INPUT .node file] [INPUT .top file] [new parameter radius] [central node] [density] [shear_visc] [bulk_visc] [shear_mod] [bulk_mod] [dielec]\n")

# Get args
inmatfname = sys.argv[1]
outmatfname = sys.argv[2]
innodefname = sys.argv[3]
intopfname = sys.argv[4]
radius = float(sys.argv[5])
node_index = int(sys.argv[6])
d = float(sys.argv[7])
sv = float(sys.argv[8])
bv = float(sys.argv[9])
sm = float(sys.argv[10])
bm = float(sys.argv[11])
di = float(sys.argv[12])

# Build objects
mat = FFEA_material.FFEA_material(inmatfname)
node = FFEA_node.FFEA_node(innodefname)
top = FFEA_topology.FFEA_topology(intopfname)

# Run sweep and set parameters
central_node = node.pos[node_index]
num_changed = 0
elindex = -1
for el in top.element:
	elindex += 1
	if np.linalg.norm(el.calc_centroid(node) - central_node) < radius:
		mat.set_params(elindex, d, sv, bv, sm, bm, di)
		num_changed += 1

print "Changed " + str(num_changed) + " elements, leaving " + str(top.num_elements - num_changed) + " elements as they were."
mat.write_to_file(outmatfname)
	
