cmake_minimum_required(VERSION 2.8.5)
project(FFEA CXX)
cmake_policy(VERSION 2.8)

# set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
set(PACKAGE_NAME "emppack")
set(VERSION "0")
set(MINOR   "9")
set(PACKAGE_VERSION "${VERSION}.${MINOR}")


option(USE_OPENMP "Build with OpenMP support" ON)


set (LIBS_TO_LINK)
#################################
#### LIBRARIES AND CFLAGS #######

# 1 #### OPENMP #####
if (USE_OPENMP)
  add_definitions(-DUSE_OPENMP)
  set(USE_PARALLEL ON)
  include(${CMAKE_ROOT}/Modules/FindOpenMP.cmake)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()
    message(FATAL_ERROR "Could not find a CXX_FLAG for OpenMP. Possible
solutions:
  -- Disable OPENMP option (-DOPENMP=0)
  -- Provide valid compiler flags for variable: OpenMP_CXX_FLAGS")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)

# 2 ### MATH ######
if (UNIX)
  set(LIBS_TO_LINK ${LIBS_TO_LINK} "m") 
endif(UNIX)

######## END OF LIBRARIES #######
#################################



#################################
###### SOURCES TO BUILD  ########
add_subdirectory("${PROJECT_SOURCE_DIR}/src")
