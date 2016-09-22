
add_test(NAME steric_n_springs COMMAND ${PROJECT_BINARY_DIR}/src/ffea cube_x.ffea)
set_tests_properties(steric_n_springs PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)

include_directories("${PROJECT_SOURCE_DIR}/include")
add_executable(test_steric checkIntersections.cpp
               ${PROJECT_SOURCE_DIR}/src/mat_vec_fns_II.cpp
               ${PROJECT_SOURCE_DIR}/src/BlobLite.cpp
               ${PROJECT_SOURCE_DIR}/src/VolumeIntersection.cpp
               ${PROJECT_SOURCE_DIR}/src/FFEA_return_codes.cpp)

add_test(NAME steric_n_springs_check COMMAND 
         ${CMAKE_CURRENT_BINARY_DIR}/test_steric 2e-9 ../cube_springs_structure/cube_x.node ../cube_springs_structure/cube_x.top 2e-9 ../cube_springs_structure/cube_x.node ../cube_springs_structure/cube_x.top cube_x_trajectory.ftj 10)
set_tests_properties(steric_n_springs_check PROPERTIES DEPENDS "steric_n_springs")




