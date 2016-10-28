cmake_minimum_required(VERSION 3.6)

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

cmake_minimum_required(VERSION 2.8)
project(triqs_interface CXX)
set(CMAKE_BUILD_TYPE Release)

include(${TRIQS_PATH}/share/triqs/cmake/TRIQSConfig.cmake)
#enable_testing()

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

#FIX ME: I want to load include dirs for HDF5 from TRIQS
find_package(HDF5 REQUIRED)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Add the compiling options (-D... ) for C++
add_definitions(${TRIQS_CXX_DEFINITIONS})

#Add include dirs
include_directories(${TRIQS_PATH}/include ${TRIQS_INCLUDE_BOOST} ${TRIQS_INCLUDE_MPI} ${TRIQS_INCLUDE_PYTHON} ${TRIQS_INCLUDE_HDF5} ${HDF5_INCLUDE_DIRS})

#Set link libraries
set(LIB_FILES ./src/solver_core.cpp)
ADD_LIBRARY(triqs_interface ${LIB_FILES})
