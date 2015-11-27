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
option(USE_FAST "Add automatic compiler flags" OFF)
mark_as_advanced(USE_FAST)
option(USE_DOUBLE "Use double precision" ON)
mark_as_advanced(USE_DOUBLE)


set (LIBS_TO_LINK)
#################################
#### LIBRARIES AND CFLAGS #######
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
set (Boost_NO_SYSTEM_PATHS ON)
set (Boost_NO_BOOST_CMAKE ON)
find_package( Boost COMPONENTS program_options )
if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  set(LIBS_TO_LINK ${LIBS_TO_LINK} ${Boost_LIBRARIES})
else (Boost_FOUND)
  message(FATAL_ERROR "Boost not found. Please provide valid paths for either:
   -- BOOST_ROOT or
   -- BOOST_INCLUDEDIR and BOOST_LIBRARYDIR")
endif (Boost_FOUND)

# 6 ### EIGEN3 #####
set(EIGEN3_HOME ${EIGEN3_HOME} CACHE PATH "[OPTIONAL] Path to your local EIGEN3 installation")
if (EIGEN3_HOME)
  find_path(EIGEN3_I_OK Sparse PATHS "${EIGEN3_HOME}/include/eigen3/Eigen"
	  "${EIGEN3_HOME}/include/Eigen" NO_DEFAULT_PATH)
  if (EIGEN3_I_OK)
    set(EIGEN3_INC "${EIGEN3_I_OK}" CACHE PATH "Path to the Eigen3 headers" FORCE)
  else (EIGEN3_I_OK)
    message(FATAL_ERROR "Eigen/Sparse not found in ${EIGEN3_HOME}/include/eigen3")
  endif (EIGEN3_I_OK)
  unset(EIGEN3_I_OK CACHE)
else(EIGEN3_HOME)
  find_path(EIGEN3_INC "Sparse" HINTS
    $ENV{EIGEN_HOME}/include/eigen3/Eigen
    $ENV{EIGEN3_HOME}/include/eigen3/Eigen
    $ENV{EIGEN_HOME}/include/Eigen
    $ENV{EIGEN3_HOME}/include/Eigen
    ${CMAKE_INCLUDE_PATH}/Eigen 
    ${CMAKE_INCLUDE_PATH}/eigen3/Eigen 
    ${CMAKE_PREFIX_PATH}/include/eigen3/Eigen
    ${CMAKE_PREFIX_PATH}/include/Eigen
    /usr/include/eigen3/Eigen
    /usr/include/Eigen
    /usr/local/include/eigen3/Eigen
    /usr/local/include/Eigen
    /opt/local/include/eigen3/Eigen
    /opt/local/include//Eigen
  )

  if(NOT EIGEN3_INC)
    message(FATAL_ERROR "Could NOT find EIGEN3: use variable EIGEN3_HOME if you have any local installation")
  endif(NOT EIGEN3_INC)
endif(EIGEN3_HOME)
set (EIGEN3_INC "${EIGEN3_INC}/..")
message(STATUS "Found Eigen3 headers: ${EIGEN3_INC}")
include_directories("${EIGEN3_INC}")


# 6 ### DOXYGEN ##### 
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file("${PROJECT_SOURCE_DIR}/Doxygen.in"
                 "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE}
                        ${PROJECT_BINARY_DIR}/Doxyfile
                        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/src"
                        COMMENT "Generating documentation for developers"
                        VERBATIM)
  install(FILES "${PROJECT_BINARY_DIR}/src/doc/html" DESTINATION
                "${CMAKE_INSTALL_PREFIX}/share/doc" OPTIONAL)
endif(DOXYGEN_FOUND)

# 7 ### DEBUG FLAGS ####
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")


# 8 ### FAST! ##########
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
add_subdirectory("${PROJECT_SOURCE_DIR}/FFEA_tools")
