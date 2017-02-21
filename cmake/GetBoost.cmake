# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

function(get_boost_toolset BOOST_TOOLSET) 
   set(BOOST_TOOLSET "UNKNOWN")
   if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
      if (APPLE)
        set(BOOST_TOOLSET "darwin" PARENT_SCOPE)
      else(APPLE)
        set(BOOST_TOOLSET "gcc" PARENT_SCOPE)
      endif(APPLE)
   endif()

   if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      set(BOOST_TOOLSET "intel" PARENT_SCOPE)
   endif()

   if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
      set(BOOST_TOOLSET "clang" PARENT_SCOPE)
   endif()

   if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Embarcadero, Borland")
      set(BOOST_TOOLSET "borland" PARENT_SCOPE)
   endif()

   if (${CMAKE_CXX_COMPILER_ID} STREQUAL "HP")
      set(BOOST_TOOLSET "hp_cxx" PARENT_SCOPE)
   endif()

   if (${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
      set(BOOST_TOOLSET "msvc" PARENT_SCOPE)
   endif()

   if (${CMAKE_CXX_COMPILER_ID} STREQUAL "SunPro")
      set(BOOST_TOOLSET "sun" PARENT_SCOPE)
   endif()
endfunction()

# 1 - Get the toolset:
get_boost_toolset(BOOST_TOOLSET)
if (BOOST_TOOLSET STREQUAL "UNKNOWN")
   message(FATAL_ERROR "I don't know how to build Boost... you'd better do that yourself, and reconfigure with USE_BOOST_EXTERNAL")
endif (BOOST_TOOLSET STREQUAL "UNKNOWN")

# 2 - Set compiler flags:
set(BOOST_COMPILER_FLAGS ${USE_CXX11_COMPILER_FLAGS})

# 3 - Set some folders:
set(BOOST_SOURCES "${PROJECT_SOURCE_DIR}/external/boost")
set(BOOST_BUILD_DIR "${PROJECT_BINARY_DIR}/external/build_boost")
set(BOOST_INSTALL_DIR "${PROJECT_BINARY_DIR}/external/install_boost")

# 4 - Configure user-config.jam
configure_file("${BOOST_SOURCES}/tools/build/user-config.jam.in"
               "${BOOST_BUILD_DIR}/build_boost/user-config.jam" @ONLY)

# --user-config="${PROJECT_BINARY_DIR}/build_boost/user-config.jam"
# 5 - Configure, Build and Install: 
include(ExternalProject)
ExternalProject_Add(ourboost
                    SOURCE_DIR ${BOOST_SOURCES}
                    INSTALL_DIR ${BOOST_INSTALL_DIR}
                    BUILD_IN_SOURCE 1
                    BUILD_ALWAYS 0
                    DEPENDERS ffea test_script_loader
                    CONFIGURE_COMMAND ""
                    BUILD_COMMAND ${BOOST_INSTALL_DIR}/bin/b2 link=static toolset=${BOOST_TOOLSET} --build-dir=${BOOST_BUILD_DIR} install --prefix=${BOOST_INSTALL_DIR} --user-config=${BOOST_BUILD_DIR}/build_boost/user-config.jam
                    INSTALL_COMMAND ""
                   )

ExternalProject_Add_Step(ourboost unix_configure DEPENDERS configure
                         WORKING_DIRECTORY ${BOOST_SOURCES}/tools/build
                         COMMAND ./bootstrap.sh --with-toolset=${BOOST_TOOLSET}
                         COMMAND ./b2 -j 4 toolset=${BOOST_TOOLSET} --build-dir=${BOOST_BUILD_DIR} install --prefix=${BOOST_INSTALL_DIR} --user-config=${BOOST_BUILD_DIR}/build_boost/user-config.jam)

set(Boost_LIBRARY_DIRS ${BOOST_INSTALL_DIR}/lib)
set(Boost_INCLUDE_DIRS ${BOOST_INSTALL_DIR}/include)
set(Boost_LIBRARIES  ${BOOST_INSTALL_DIR}/lib/libboost_program_options.a;${BOOST_INSTALL_DIR}/lib/libboost_filesystem.a;${BOOST_INSTALL_DIR}/lib/libboost_system.a)


