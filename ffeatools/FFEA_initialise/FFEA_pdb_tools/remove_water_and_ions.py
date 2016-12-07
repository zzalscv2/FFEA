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
import FFEA_pdb

if len(sys.argv) != 3:
	sys.exit("Usage: python convert_dynein.py [INPUT .pdb fname] [OUTPUT .pdb fname]")

# Get args
infname = sys.argv[1]
outfname = sys.argv[2]

# Do stuff
pdb = FFEA_pdb.FFEA_pdb(infname)
pdb.reorder_atoms()
pdb.reorder_residues()
pdb.write_to_file(outfname)
