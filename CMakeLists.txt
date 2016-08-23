cmake_minimum_required(VERSION 2.8.5)
project(FFEA CXX C)
cmake_policy(VERSION 2.8)

# set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
set(PACKAGE_NAME "emppack")
set(VERSION "0")
set(MINOR   "9")
set(PACKAGE_VERSION "${VERSION}.${MINOR}")


option(USE_OPENMP "Build with OpenMP support" ON)
set(USE_OMP_MODE 1 CACHE STRING "0: disabled; 1: threads within blob; 2: threads per blob")
option(USE_MPI "Use MPI parallelisation: EXPERIMENTAL!" OFF)
mark_as_advanced(USE_MPI)
option(USE_FAST "Add automatic compiler flags" ON)
mark_as_advanced(USE_FAST)
option(USE_DOUBLE "Use double precision" ON)
mark_as_advanced(USE_DOUBLE)


set (LIBS_TO_LINK)
#################################
#### LIBRARIES AND CFLAGS #######

# benchmark ### MPI #########
if (USE_MPI)
  add_definitions(-DUSE_MPI)

  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
  include_directories(${MPI_INCLUDE_PATH})
  set(LIBS_TO_LINK ${LIBS_TO_LINK} ${MPI_LIBRARIES}) 
endif(USE_MPI)


# 0 #### PRECISION ###
if (USE_DOUBLE)
  add_definitions(-DUSE_DOUBLE)
endif(USE_DOUBLE)

# 1 #### OPENMP #####
if (USE_OPENMP)
  add_definitions(-DUSE_OPENMP)
endif(USE_OPENMP)
include(${CMAKE_ROOT}/Modules/FindOpenMP.cmake)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_FOUND)
    if (USE_OPENMP)
      message(FATAL_ERROR "Could not find a CXX_FLAG for OpenMP. Possible
solutions:
  -- Disable OPENMP option (-DUSE_OPENMP=0)
  -- Provide valid compiler flags for variable: OpenMP_CXX_FLAGS")
    else (USE_OPENMP)
      message(FATAL_ERROR "Could not find a CXX_FLAG for OpenMP.  Possible
solutions:
  -- Provide valid compiler flags for variable: OpenMP_CXX_FLAGS
  -- Submit a request to the developers suggesting to change omp_get_wtime to something else")
    endif (USE_OPENMP)
endif(OPENMP_FOUND)

# 2 #### PARALLEL MODE #########
if (USE_OMP_MODE EQUAL "0")
  if (USE_OPENMP) 
    message(FATAL_ERROR "OMP was enabled, but USE_OMP_MODE set to 0. Possible solutions: 
    -- Disable OPENMP option (-DUSE_OPENMP=0) or
    -- Set USE_OMP_MODE either to 1 or 2")
  endif (USE_OPENMP)
else ()
  if (USE_OMP_MODE EQUAL "1")
   add_definitions(-DFFEA_PARALLEL_WITHIN_BLOB)
  else ()
    if (USE_OMP_MODE EQUAL "2")
      add_definitions(-DFFEA_PARALLEL_PER_BLOB)
    else ()
      message(FATAL_ERROR "invalid value for USE_OMP_MODE")
    endif ()
  endif()
endif()

# 3 #### USE C++11 ############ 
# g++: -std=c++0x / -std=c++11
# icpc: -std=c++0x / -std=c++11 on LINUX, OS X and WINDOWS 
# pg++: 13.2 crashes with flag -std=c++0x, besides it does not support
#         the latest standard ("cannot open source file "future"). 
# MSVC: does not need any flag.
include(cmake/CheckForCXX11Features.cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11_COMPILER_FLAGS}")


# 4 ### MATH ######
if (UNIX)
  set(LIBS_TO_LINK ${LIBS_TO_LINK} "m") 
endif(UNIX)

# 5 ### BOOST #####
# set (Boost_NO_SYSTEM_PATHS ON)
set (Boost_NO_BOOST_CMAKE ON)
# Check version:
find_package( Boost 1.54.0 REQUIRED COMPONENTS program_options filesystem system )
if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -L${Boost_LIBRARY_DIRS}")
  set(LIBS_TO_LINK ${LIBS_TO_LINK} ${Boost_LIBRARIES})
else (Boost_FOUND)
  message(FATAL_ERROR "You could provide valid paths through either:
   -- BOOST_ROOT or
   -- BOOST_INCLUDEDIR and BOOST_LIBRARYDIR")
endif (Boost_FOUND)

# 6 ### EIGEN3 #####
set(EIGEN3_HOME $ENV{EIGEN3_HOME} CACHE PATH "[OPTIONAL] Path to your local EIGEN3 installation")
# Sometimes "FindEigen3.cmake" is not distributed with CMake, so we're shipping it:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
find_package(Eigen3 3.2.0 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

# 6 ### PYTHON ######
find_package(PythonInterp REQUIRED)
# python will be used to run some tests.


# 7 ### DOXYGEN ##### 
find_package(Doxygen 1.8)
if(DOXYGEN_FOUND)
  configure_file("${PROJECT_SOURCE_DIR}/Doxygen.in"
                 "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE}
                        ${PROJECT_BINARY_DIR}/Doxyfile
                        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/src"
                        COMMENT "Generating documentation for developers"
                        VERBATIM)
  install(DIRECTORY "${PROJECT_BINARY_DIR}/src/doc/html" DESTINATION
                "${CMAKE_INSTALL_PREFIX}/share/ffea/doc" OPTIONAL)
endif(DOXYGEN_FOUND)

# 8 ### DEBUG FLAGS ####
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")


# 9 ### FAST! ##########
if (USE_FAST)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -ffast-math -fipa-pta")
    else(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -flto -ffast-math -fipa-pta -march=native")
    endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ipo -no-prec-div -fp-model fast=2 -xHost")
  endif()
  message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")
endif(USE_FAST)

######## END OF LIBRARIES #######
#################################



#############################################
###### SOURCES TO BUILD AND INSTALL  ########
add_subdirectory("${PROJECT_SOURCE_DIR}/src")


#######################################
######## INSTALL THE FFEA TOOLS  ######
# 1 - Set up the install python path
include(${CMAKE_ROOT}/Modules/FindPythonInterp.cmake)
set(PYTHONSTUFF "${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")

# 2 - and go through the folders:
# PENDING # Unify or rename make_structure_map.cpp Matrices.hpp and Vectors.hpp in
#               FFEA_mapping_tools and FFEA_node_tools
# PENDING # FFEA_initialise/PDB_tools # deprecate FFEA_convert_traj_to_pdb.c (?)
# EVEN WORSE: make_structure_map2 in FFEA_initialise/FFEA_node_tools
add_subdirectory("${PROJECT_SOURCE_DIR}/ffeatools")



#####################################################
############ # # # ENABLE TESTS # # # ###############
enable_testing()
# set_tests_properties(ENVIRONMENT OMP_NUM_THREADS=1)
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/consistency")
# add_subdirectory("${PROJECT_SOURCE_DIR}/tests/consistency/volume")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/physics")
