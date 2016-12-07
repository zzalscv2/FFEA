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

add_executable(FFEA_extract_positions FFEA_extract_positions.c)
add_executable(FFEA_extract_specific_nodes FFEA_extract_specific_nodes.c)

install(TARGETS FFEA_extract_positions FFEA_extract_specific_nodes
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis/FFEA_traj_tools")

# install(TARGETS FFEA_extract_positions FFEA_extract_specific_nodes
#         DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

install(FILES __init__.py
        convert_old_traj_to_new.py Dx_convert_from_xyz_position_trajectory.py
        Dx_smooth_convert_from_xyz_position_trajectory.py FFEA_convert_old_traj_to_new.py
        FFEA_map_trajectory_to_PDB.py FFEA_thin_trajectory.py FFEA_traj_to_nodes.py
        FFEA_traj_to_PDB_traj.py FFEA_convert_traj_to_pdb.py FFEA_trim_trajectory.py FFEA_split_trajectory.py
        FFEA_get_snapshots_in_nodes.py FFEA_strip_equilibration.py FFEA_get_num_frames.py PDB_convert_to_FFEA_trajectory.py
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis/FFEA_traj_tools")

