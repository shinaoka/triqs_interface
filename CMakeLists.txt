# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# start configuration
cmake_minimum_required(VERSION 2.8)
project(alps_cthyb CXX)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Add the compiling options (-D... ) for C++
message(STATUS "TRIQS : Adding compilation flags detected by the library (C++11/14, libc++, etc...) ")
add_definitions(${TRIQS_CXX_DEFINITIONS})

# Options for the compilation ...
#option(Tests "Enable Tests" ON)

include_directories(${CMAKE_BINARY_DIR}/c++ c++)

# Compile C++ code
add_subdirectory(c++)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)
 if (${Tests})
 add_subdirectory(test)
 endif()
 option(BUILD_DOC "Build documentation" OFF)
 if(${BUILD_DOC})
   if(NOT ${TRIQS_WITH_DOCUMENTATION})
     message(FATAL_ERROR "Error: TRIQS library has not been compiled with its documentation")
   endif()
   add_subdirectory(doc)
 endif(${BUILD_DOC})
endif()