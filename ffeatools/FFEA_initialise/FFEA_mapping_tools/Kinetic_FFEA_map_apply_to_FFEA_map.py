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
import numpy as np
import FFEA_kinetic_map

if len(sys.argv) != 5:
	sys.exit("Usage: python " + os.path.basename(os.path.abspath(sys.argv[0])) + " [INPUT MAP A] [INPUT MAP B] [OUTPUT AB] [OUTPUT BA]")

# Get args
mapa_fname = sys.argv[1]
mapb_fname = sys.argv[2]
outmapab_fname = sys.argv[3]
outmapba_fname = sys.argv[4]

# Get maps
mapa = FFEA_kinetic_map.FFEA_kinetic_map(mapa_fname)
mapb = FFEA_kinetic_map.FFEA_kinetic_map(mapb_fname)

# Apply maps
mapab = mapa.apply_to_map(mapb)
mapba = mapb.apply_to_map(mapa)

# Write maps
mapab.write_to_file(outmapab_fname)
mapba.write_to_file(outmapba_fname)
