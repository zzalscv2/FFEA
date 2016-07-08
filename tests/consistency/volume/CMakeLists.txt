include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(test_volume testVolume.cpp
               ${PROJECT_SOURCE_DIR}/src/mat_vec_fns_II.cpp
               ${PROJECT_SOURCE_DIR}/src/VolumeIntersection.cpp)

add_test(test_volume ${CMAKE_CURRENT_BINARY_DIR}/test_volume)
