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
import FFEA_surf
if len(sys.argv) != 4:
	sys.exit("Usage: python scale_netgen_surf.py [INPUT .surf] [OUTPUT .surf] [scale]")

# Get args
insurffname = sys.argv[1]
outsurffname = sys.argv[2]
scale = float(sys.argv[3])

# Make and scale a surf
surf = FFEA_surf.FFEA_surf(insurffname)
surf.scale(scale)

# Write to out_fname
surf.write_to_netgen_surf(outsurffname)


