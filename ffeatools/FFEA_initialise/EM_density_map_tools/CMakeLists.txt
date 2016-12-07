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

add_executable(add_maps add_maps.cpp)
add_executable(convert_pdb_to_emdb_map convert_pdb_to_emdb_map.cpp)

install (TARGETS add_maps convert_pdb_to_emdb_map RUNTIME 
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/EM_density_map_tools")

install (FILES __init__.py convert_pdb_to_emdb_map.py
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/EM_density_map_tools")

