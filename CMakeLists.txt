cmake_minimum_required(VERSION 2.8.5)
project(FFEA CXX C)
cmake_policy(VERSION 2.8)

# set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
set(PACKAGE_NAME "ffea")
set(VERSION "1")
set(MINOR   "0")
set(PACKAGE_VERSION "${VERSION}.${MINOR}")


option(USE_OPENMP "Build with OpenMP support" ON)
set(USE_OMP_MODE 1 CACHE STRING "0: disabled; 1: threads within blob; 2: threads per blob")
option(USE_MPI "Use MPI parallelisation: EXPERIMENTAL!" OFF)
mark_as_advanced(USE_MPI)
option(USE_FAST "Add automatic compiler flags" ON)
mark_as_advanced(USE_FAST)
option(USE_PACK "Add genereic automatic compiler flags" OFF)
mark_as_advanced(USE_PACK)
option(USE_PACK_DIST "accepted values are RPM and DEB" "RPM")
mark_as_advanced(USE_PACK_DIST)
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
set (Boost_USE_STATIC_LIBS ON)
# Check version:
find_package( Boost 1.54.0 REQUIRED COMPONENTS program_options filesystem system  )
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
# 6.1 # set up EIGEN3_HOME
if (EIGEN3_HOME)  
  set(EIGEN3_HOME ${EIGEN3_HOME} CACHE PATH "[OPTIONAL] Path to your local EIGEN3 installation")
else(EIGEN3_HOME)   
  set(EIGEN3_HOME $ENV{EIGEN3_HOME} CACHE PATH "[OPTIONAL] Path to your local EIGEN3 installation")
endif(EIGEN3_HOME)
# 6.2 # we guarantee that if $EIGEN3_HOME is given, it will be take 
      # before than any of the other possible installations available.
set(CMAKE_PREFIX_PATH ${EIGEN3_HOME} $ENV{EIGEN_HOME} ${CMAKE_PREFIX_PATH})
# 6.3 # Sometimes "FindEigen3.cmake" is not distributed with CMake, so we're shipping it:
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
# 6.4 # Go for it:
find_package(Eigen3 3.2.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})



# 7 ### PYTHON ######
find_package(PythonInterp REQUIRED)
# python will be used to run some tests.


# 8 ### DOXYGEN ##### 
find_package(Doxygen 1.8)
if(DOXYGEN_FOUND)
  configure_file("${PROJECT_SOURCE_DIR}/Doxygen.in"
                 "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
  configure_file("${PROJECT_SOURCE_DIR}/DoxygenPyM.in"
                 "${PROJECT_BINARY_DIR}/DoxyfilePyM" @ONLY)
  add_custom_target(docpy ALL ${DOXYGEN_EXECUTABLE}
                        ${PROJECT_BINARY_DIR}/DoxyfilePyM
                        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/src"
                        COMMENT "Generating documentation for developers"
                        VERBATIM)
  add_custom_target(doc ALL ${DOXYGEN_EXECUTABLE}
                        ${PROJECT_BINARY_DIR}/Doxyfile
                        DEPENDS docpy
                        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/src"
                        COMMENT "Generating documentation for developers"
                        VERBATIM)
  install(DIRECTORY "${PROJECT_BINARY_DIR}/src/doc/html" DESTINATION
                share/ffea/doc OPTIONAL)
  install(DIRECTORY "${PROJECT_BINARY_DIR}/src/ffeamodules/html" DESTINATION
                share/ffea/ffeamodules OPTIONAL)
endif(DOXYGEN_FOUND)

# 9 ### DEBUG FLAGS ####
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")


# 10 ### FAST! ##########
if (NOT USE_PACK)
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
else (NOT USE_PACK)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -flto -ffast-math -fipa-pta -mtune=generic")
  else()
      message(FATAL_ERROR "No packaging without GNU compilers")
  endif()
  if (USE_FAST)
      message(STATUS "USE_FAST changed to USE_PACK!")
  endif(USE_FAST)
  message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")
endif(NOT USE_PACK)

######## END OF LIBRARIES #######
#################################



#############################################
###### SOURCES TO BUILD AND INSTALL  ########
add_subdirectory("${PROJECT_SOURCE_DIR}/src")
#############################################

#######################################
######## INSTALL THE FFEA TOOLS  ######
# 1 - Set up the install python path
include(${CMAKE_ROOT}/Modules/FindPythonInterp.cmake)
option(USE_PACK_DIST "accepted values are RPM and DEB" "RPM")
set(PYTHONSTUFF "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/ffeatools")
# Debian and derivatives use dist-packages instead of site-packages. 
if (USE_PACK_DIST STREQUAL DEB)
  set(PYTHONSTUFF "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/dist-packages/ffeatools")
endif(USE_PACK_DIST STREQUAL DEB)

# 2 - and go through the folders:
add_subdirectory("${PROJECT_SOURCE_DIR}/ffeatools")
#######################################



#####################################################
############ # # # ENABLE TESTS # # # ###############
enable_testing()
# set_tests_properties(ENVIRONMENT OMP_NUM_THREADS=1)
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/consistency")
# add_subdirectory("${PROJECT_SOURCE_DIR}/tests/consistency/volume")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/physics")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/io")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/ffeatools")
#####################################################

###############################
##### INSTALL THE LICENSE #####
install(FILES COPYING DESTINATION share/ffea COMPONENT data)

#####################################################
########### # # # CPACK STUFF # # # #################
# FFEA links to Boost statically,
#  so if packagising for different platforms, we only need their gcc version: 
#  gcc 4.8 - Native - Ubuntu 14.04, Mint 17, Centos 7
#          - Works on - Ubuntu 16.04
#  gcc 5.4 - needs Eigen 3.2.10 - 
#          - Ubuntu 16.04 
#  gcc 4.4 - CentOS 6
set(CPACK_PACKAGE_NAME ${PACKAGE_NAME})
set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION})
set(CPACK_PACKAGE_ARCHITECTURE "x86_64")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_${CPACK_PACKAGE_ARCHITECTURE}")
set(CPACK_PACKAGE_VENDOR "University of Leeds, UK.")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "FFEA - a toolkit for molecular simulation using Fluctuating Finite Element Method")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_CONTACT "contact-us@leeds.ac.uk")

# # # CPACK: DEB or RPM # # # 
message(STATUS "Got DIST ${USE_PACK_DIST}")
if (USE_PACK_DIST STREQUAL "DEB") 
  set(CPACK_GENERATOR "TBZ2" "DEB")
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_amd64")
elseif (USE_PACK_DIST STREQUAL "RPM")
  set(CPACK_GENERATOR "TBZ2" "RPM")
endif()

# # # CPACK: DEB attributes # # # 
# Choose: 
# sudo dpkg -i ffea_1.0_amd64.deb
# sudo apt-get install -f 
# or 
# sudo apt-get install python2.7 python-numpy python-matplotlib libgomp1
# sudo dpkg -i ffea_1.0_amd64.deb
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
set(CPACK_DEBIAN_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libgomp1 (>=5.4.0), python2.7, python-numpy, python-matplotlib")
set(CPACK_DEBIAN_PACKAGE_RECOMMENDS "pymol (>=1.8)")

# # # CPACK: RPM attributes # # # 
# sudo yum install numpy python-matplotlib
# sudo rpm -Uvh ffea_1.0_x86_64.rpm
set(CPACK_RPM_PACKAGE_LICENSE "GPLv3")
set(CPACK_RPM_PACKAGE_REQUIRES "numpy, python-matplotlib")
# pymol 1.6 in CentOS 6 is too old for our plugin and CentOS 7 does not provide pymol.


#must come after all cpack settings!
if (USE_PACK)
  include(CPack)
endif(USE_PACK)

