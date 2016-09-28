FIND_PACKAGE(Boost REQUIRED)

set (TESTSCRIPT "${PROJECT_BINARY_DIR}/tests/io/script")

include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(test_script_loader test_script_loader.cpp
               ${PROJECT_SOURCE_DIR}/src/FFEA_input_reader.cpp
               ${PROJECT_SOURCE_DIR}/src/SimulationParams.cpp
               ${PROJECT_SOURCE_DIR}/src/FFEA_return_codes.cpp)

file (COPY 2r5u_8ang.ffea DESTINATION ${TESTSCRIPT})

target_link_libraries(test_script_loader ${LIBS_TO_LINK})

add_test(test_script_loader ${CMAKE_CURRENT_BINARY_DIR}/test_script_loader)

