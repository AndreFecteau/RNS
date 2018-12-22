# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required (VERSION 2.4 FATAL_ERROR)

# Set the project name
project (RNS_Solver CXX)

# project version
set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

# detect operating system and host processor
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")

# Set a default build type for single-configuration
# CMake generators if no build type is set.
if(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

# Select flags.
SET(CMAKE_CXX_FLAGS "-Wall -Wextra")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG -march=native")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

# try conditional compilation
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

# check results and add flag
if(COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ldl -pthread -std=c++17")
    set(CMAKE_CXX_STANDARD 17)
    message(STATUS "Using the compiler C++17 ${CMAKE_CXX_COMPILER}")
elseif(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ldl -pthread -std=c++14")
    set(CMAKE_CXX_STANDARD 14)
    message(STATUS "Using the compiler C++14 ${CMAKE_CXX_COMPILER}")
elseif(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ldl -pthread -std=c++11")
    set(CMAKE_CXX_STANDARD 11)
    message(STATUS "Using the compiler C++11 ${CMAKE_CXX_COMPILER}")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 || C++14 || C++17 support. Please use a different C++ compiler.")
endif()

option(USE_OPENMP "Use OpenMP in build" ON)

# Check user include dirmake
  include_directories (SYSTEM "~/include")
# Standard path for libraries
  include_directories (SYSTEM "/usr/include/eigen3")
  include_directories (SYSTEM "/usr/include/cereal")

if (USE_OPENMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopenmp")
  message(STATUS "Linked OpenMP compiler flag")
endif (USE_OPENMP)

add_executable (RNS_Solver RNS_Main.cpp)

message(STATUS "Generated build with CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
