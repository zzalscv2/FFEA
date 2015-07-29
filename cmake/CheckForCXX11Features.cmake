
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
           "${CMAKE_CURRENT_BINARY_DIR}/cxx11_stoi}"
           "${PROJECT_SOURCE_DIR}/cmake/tests/stoi.cpp"
           COMPILE_DEFINITIONS "${CXX11_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}"
           OUTPUT_VARIABLE OUTPUT)

if (RESULT_VAR)
  message(STATUS "Checking C++11 support for stoi - Success")
else ()
  message(FATAL_ERROR "Checking C++11 support for stoi - Failed")
endif ()
