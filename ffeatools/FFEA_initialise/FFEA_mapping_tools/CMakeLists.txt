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

add_executable(make_structure_map make_structure_map.cpp)

install (TARGETS make_structure_map 
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_mapping_tools")

# install (TARGETS make_structure_map 
#          DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")


install(FILES FFEA_map_PDB_to_FFEA.py FFEA_convert_kinetic_map_to_sparse.py FFEA_make_structures_overlap.py FFEA_generate_kinetic_maps.py Kinetic_FFEA_map_apply_to_FFEA_node.py Kinetic_FFEA_map_apply_to_FFEA_traj.py Kinetic_FFEA_map_apply_to_FFEA_map.py Kinetic_FFEA_map_apply_to_PDB.py
         DESTINATION "${PYTHONSTUFF}/FFEA_initialise/FFEA_mapping_tools")
