# copy every folder here to PROJECT_BINARY_DIR to run the tests. 
#   consider configuring the number of steps. 
set (TPHYSICS "${PROJECT_BINARY_DIR}/tests/physics")
file (MAKE_DIRECTORY ${TPHYSICS})
file (COPY cuboid_flat cuboid_flat_structure sphere_63_120_structure sphere_63_120_two
                  sphere_63_120_two_vdw sphere_63_120_cuboid_flat sphere_63_120_mass 
                  sphere_63_120_nomass sphere_63_120_nomass_restart
                  sphere_63_120_nomass_two_vdw-preComp 
                  cube_springs_structure steric_cubes_w_springs
                  DESTINATION ${TPHYSICS})

# and run the tests. 
add_subdirectory(sphere_63_120_mass)
add_subdirectory(sphere_63_120_nomass)
add_subdirectory(sphere_63_120_nomass_restart)
add_subdirectory(steric_cubes_w_springs)
