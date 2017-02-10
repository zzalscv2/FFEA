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


include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" _HAS_CXX11_FLAG)
if (NOT _HAS_CXX11_FLAG)
    check_cxx_compiler_flag("-std=c++0x" _HAS_CXX0X_FLAG)
endif ()

if (_HAS_CXX11_FLAG)
    set(CXX11_COMPILER_FLAGS "-std=c++11")
elseif (_HAS_CXX0X_FLAG)
    set(CXX11_COMPILER_FLAGS "-std=c++0x")
endif ()

##############################
# CHECK 1 #### Check for stoi: 
message(STATUS "Checking C++11 support for stoi")
try_compile(RESULT_VAR 
           "${CMAKE_CURRENT_BINARY_DIR}/cxx11_stoi"
           "${PROJECT_SOURCE_DIR}/cmake/tests_compilation/stoi.cpp"
           COMPILE_DEFINITIONS "${CXX11_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}"
           OUTPUT_VARIABLE OUTPUT)

if (RESULT_VAR)
  message(STATUS "Checking C++11 support for stoi - Success")
else ()
  message(FATAL_ERROR "Checking C++11 support for stoi - Failed")
endif ()


##############################
# CHECK 2 #### Check for fstream open(string): 
message(STATUS "Checking C++11 support for fstream-string")
try_compile(RESULT_VAR 
           "${CMAKE_CURRENT_BINARY_DIR}/cxx11_fstream-string"
           "${PROJECT_SOURCE_DIR}/cmake/tests_compilation/fstream-string.cpp"
           COMPILE_DEFINITIONS "${CXX11_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}"
           OUTPUT_VARIABLE OUTPUT)

if (RESULT_VAR)
  message(STATUS "Checking C++11 support for fstream-string - Success")
else ()
  message(FATAL_ERROR "Checking C++11 support for fstream-string - Failed")
endif ()


##############################
# CHECK 3 #### Check for future threads: 
message(STATUS "Checking C++11 support for future threads")
set(CHECK_REQUIRED_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
try_compile(RESULT_VAR 
            "${CMAKE_CURRENT_BINARY_DIR}/cxx11_future"
            "${PROJECT_SOURCE_DIR}/cmake/tests_compilation/cxx11-test-future.cpp"
            COMPILE_DEFINITIONS "${CXX11_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}"
            CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${CHECK_REQUIRED_LIBRARIES}"
            OUTPUT_VARIABLE OUTPUT)

if (RESULT_VAR)
  message(STATUS "Checking C++11 support for future threads - Success")
  set(HAS_CXX11_FUTURE TRUE)
else ()
  message(WARNING "Checking C++11 support for future threads - Failed")
  set(HAS_CXX11_FUTURE FALSE)
endif ()
