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


add_subdirectory("Euler_characteristic_analysis")
add_subdirectory("FFEA_pyPca")
add_subdirectory("FFEA_traj_tools")

install (FILES FFEA_thin_system.py FFEA_build_analysis_environment.py __init__.py
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis")

install(DIRECTORY FFEA_meas_tools
        DESTINATION "${PYTHONSTUFF}/FFEA_analysis")

