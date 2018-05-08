# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required (VERSION 3.2 FATAL_ERROR)

# Set the project name
project (RNS_Solver CXX)

# project version
set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

# detect operating system and host processor
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")

# example how to set c++ compiler flags for GNU
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
else()
  message(FATAL_ERROR "Not GNU Compiler, not sure what will happen.")
endif()

# try conditional compilation
include(CheckCXXCompilerFlag)
# CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)

# check results and add flag
# if(COMPILER_SUPPORTS_CXX17)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
#     set(CMAKE_CXX_STANDARD 17)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    set(CMAKE_CXX_STANDARD 14)
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()

# find a Eigen install
find_package (Eigen3 3.2.8)
if(Eigen3_FOUND)
 message ("Eigen3 found")
 set(TARGET_LINK_LIBRARIES_PATH "Eigen3::Eigen")
else()
  include_directories ( "/usr/include/eigen3" )
  message ("Cannot find Eigen3; using default path /usr/include/eigen3")
endif()

# find a Cereal install
find_package (Cereal QUIET)
if(Cereal_FOUND)
 message ("Cereal found")
else()
  include_directories ( "/usr/include/cereal" )
  message ("Cannot find Cereal; using default path /usr/include/cereal")
endif()

# find a OpenMP install
find_package (OpenMP)
if(OpenMP_FOUND)
 message ("OpenMP found")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
 message (FATAL_ERROR "Cannot find OpenMP")
endif()

add_executable (RNS_Solver RNS_Main.cpp)

add_executable (RNS_Solver_restart RNS_Load.cpp)

# Link Eigen librairies
target_link_libraries (RNS_Solver ${TARGET_LINK_LIBRARIES_PATH})
target_link_libraries (RNS_Solver_restart ${TARGET_LINK_LIBRARIES_PATH})
# add_executable (RNS_Solver_test Block_Triagonal_Matrix_Inverse_Test.cpp)
# target_link_libraries (RNS_Solver_test Eigen3::Eigen)
