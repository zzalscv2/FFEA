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

import FFEA_pin, FFEA_node
import sys

if len(sys.argv) != 5:
	sys.exit("Usage: python FFEA_pin_radially.py [INPUT .node] [OUTPUT .pin] [Node to pin about] [Radius]")

# Get args
node = FFEA_node.FFEA_node(sys.argv[1])
pin = FFEA_pin.FFEA_pin(sys.argv[2])
node_index = int(sys.argv[3])
radius = float(sys.argv[4])
pin.pin_radially(node, node_index, radius)
pin.write_to_file(sys.argv[2])
