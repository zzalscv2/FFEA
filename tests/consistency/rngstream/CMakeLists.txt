include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(test_rngStream testRngStream.cpp ${PROJECT_SOURCE_DIR}/src/RngStream.cpp)

add_test(test_rngStream ${CMAKE_CURRENT_BINARY_DIR}/test_rngStream)
