cmake_minimum_required(VERSION 2.8.5)
project(FFEA CXX)
cmake_policy(VERSION 2.8)

# set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
set(PACKAGE_NAME "emppack")
set(VERSION "0")
set(MINOR   "9")
set(PACKAGE_VERSION "${VERSION}.${MINOR}")


option(USE_OPENMP "Build with OpenMP support" ON)
set(USE_OMP_MODE 1 CACHE STRING "0: disabled; 1: threads within blob; 2: threads per blob")


set (LIBS_TO_LINK)
#################################
#### LIBRARIES AND CFLAGS #######

# 1 #### OPENMP #####
if (USE_OPENMP)
  # add_definitions(-DUSE_OPENMP)
  set(USE_PARALLEL ON)
  include(${CMAKE_ROOT}/Modules/FindOpenMP.cmake)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()
    message(FATAL_ERROR "Could not find a CXX_FLAG for OpenMP. Possible
solutions:
  -- Disable OPENMP option (-DUSE_OPENMP=0)
  -- Provide valid compiler flags for variable: OpenMP_CXX_FLAGS")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)

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

# 2 ### MATH ######
if (UNIX)
  set(LIBS_TO_LINK ${LIBS_TO_LINK} "m") 
endif(UNIX)

######## END OF LIBRARIES #######
#################################



#################################
###### SOURCES TO BUILD  ########
add_subdirectory("${PROJECT_SOURCE_DIR}/src")
