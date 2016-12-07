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
import FFEA_node, FFEA_pdb

if (len(sys.argv) != 4):
	sys.exit("Usage: python FFEA_map_PDB_to_FFEA.py [INPUT .node] [INPUT .pdb] [INPUT PDB-Node scale]")
	
# Get args
innode = sys.argv[1]
inpdb = sys.argv[2]
scale = float(sys.argv[3])

# Build structures
node = FFEA_node.FFEA_node(innode)
pdb = FFEA_pdb.FFEA_pdb(inpdb)

# Assuming these are not deformed (volume overlap maximisation may be possible, but this is easier for starters)
# Move both to origin
node.set_pos([0.0,0.0,0.0])
pdb.blob[0].frame[0].set_pos([0.0,0.0,0.0])

basepdb = "temp.pdb"
targetnode = "temp.node"
mapfname = "pdbtonode.map"
sparsemapfname = "pdbtonode_sparse.map"
node.write_to_file(targetnode)
pdb.write_to_file(basepdb)

# Make the map
scriptdir = os.path.dirname(os.path.abspath(sys.argv[0]))
os.system(scriptdir + "/make_structure_map -i %s -o %s -m %s -s %f" % (basepdb, targetnode, mapfname, scale))

# Make them sparse
os.system("python " + scriptdir + "/FFEA_convert_kinetic_map_to_sparse.py " + mapfname + " " + sparsemapfname)

os.system("rm %s %s" % (basepdb, targetnode))