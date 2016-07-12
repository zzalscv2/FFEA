# copy every folder here to PROJECT_BINARY_DIR to run the tests. 
#   consider configuring the number of steps. 
set (PHYSICS "${PROJECT_BINARY_DIR}/tests/physics")
file (MAKE_DIRECTORY ${PHYSICS})
file (COPY cuboid_flat cuboid_flat_structure sphere_63_120_structure sphere_63_120_two
                  sphere_63_120_two_vdw sphere_63_120_cuboid_flat sphere_63_120_mass 
                  sphere_63_120_nomass DESTINATION ${PHYSICS})
# and run the tests. 
add_test(NAME sphere_two_vdw COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_two_vdw/sphere_63_120_two_vdw.ffea)
set_tests_properties(sphere_two_vdw PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
add_test(NAME sphere_two_vdw_check COMMAND ${CMAKE_COMMAND} -E compare_files 
              sphere_63_120_two_vdw/sphere_63_120_two_vdw_trajectory.out
              ${PROJECT_SOURCE_DIR}/tests/physics/sphere_63_120_two_vdw/sphere_63_120_two_vdw_trajectory.out)



