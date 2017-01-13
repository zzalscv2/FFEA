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

# copy every folder here to PROJECT_BINARY_DIR to run the tests. 
#   consider configuring the number of steps. 
set (TPHYSICS "${PROJECT_BINARY_DIR}/tests/physics")
file (MAKE_DIRECTORY ${TPHYSICS})
file (COPY cuboid_flat cuboid_flat_structure sphere_63_120_structure sphere_63_120_two
                  sphere_63_120_two_vdw sphere_63_120_cuboid_flat sphere_63_120_mass 
                  sphere_63_120_nomass sphere_63_120_nomass_restart
                  sphere_63_120_nomass_two_vdw-preComp 
                  cube_springs_structure steric_cubes_w_springs
		  sphere_diffusion
                  cyl_160_fine EI_cyl_160_fine E_cyl_160_fine
                  DESTINATION ${TPHYSICS})

# and run the tests. 
add_subdirectory(sphere_63_120_mass)
add_subdirectory(sphere_63_120_nomass)
add_subdirectory(sphere_63_120_nomass_restart)
add_subdirectory(sphere_63_120_nomass_two_vdw-preComp)
add_subdirectory(steric_cubes_w_springs)
add_subdirectory(sphere_diffusion)
add_subdirectory(E_cyl_160_fine)
add_subdirectory(EI_cyl_160_fine)
