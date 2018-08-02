# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# start configuration
cmake_minimum_required(VERSION 3.1)
project(alps_cthyb CXX)
if (CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()
enable_testing()

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)


# Add the compiling options (-D... ) for C++
#message(STATUS "TRIQS : Adding compilation flags detected by the library
#(C++11/14, libc++, etc...) ")
message(STATUS "Adding definitions from TRIQS: ${TRIQS_CXX_DEFINITIONS}")
message(STATUS "Please confirm that the same C++ standard was used to build
ALPSCOre libraries")
add_definitions(${TRIQS_CXX_DEFINITIONS})

# Load ALPSCore. If ALPSCore is not found
# automatically, specify its location using:
# export ALPSCore_DIR=/location/to/ALPSCORE/
find_package(ALPSCore REQUIRED)
find_package(ALPSCoreCTHYB REQUIRED)

# Check if 
if (NOT ALPSCore_CXX_STD MATCHES "^[cC][+][+]14|custom$")
    message(FATAL_ERROR "ALPSCore libraries was built with ${ALPSCore_CXX_STD}. Please rebuild ALPSCore libraries with C++14.")
endif ()

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Add the compiling options (-D... ) for C++
if (CMAKE_BUILD_TYPE STREQUAL "Release")
  add_definitions(${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -DBOOST_DISABLE_ASSERTS)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS})
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  add_definitions(${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${CMAKE_CXX_FLAGS})
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
  add_definitions(${CMAKE_CXX_FLAGS_MINSIZEREL} ${CMAKE_CXX_FLAGS})
endif()

# Options for the compilation ...

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
