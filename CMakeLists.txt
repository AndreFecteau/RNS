# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required (VERSION 3.5)


# Set the project name
project (RNS_Solver)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
  "MinSizeRel" "RelWithDebInfo")
endif()

# try conditional compilation
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)

# check results and add flag
if(COMPILER_SUPPORTS_CXX17)#
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    set(CMAKE_CXX_STANDARD 17)
elseif(COMPILER_SUPPORTS_CXX14)#
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    set(CMAKE_CXX_STANDARD 14)
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++17 or C++14 support. Please use a different C++ compiler.")
endif()

# find a Eigen install
find_package (Eigen3 3.3.4 REQUIRED NO_MODULE)

# find a Cereal install
find_package (Cereal REQUIRED)

# find a OpenMP install
find_package (OpenMP)

# Set a default C++ compile flag
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# Add OpenMP Flags if found
#if(OpenMP_FOUND)
#  message ("OpenMP found")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#else()
#  message (FATAL_ERROR "Cannot find OpenMP")
#endif()

if(Cereal_FOUND)
  message ("Cereal found")
else()
  message (FATAL_ERROR "Cannot find Cereal")
endif()

if(Eigen3_FOUND)
  message ("Eigen3 found")
else()
  message (FATAL_ERROR "Cannot find Eigen3")
endif()

# Add an executable
add_executable (RNS_Solver Setup.cpp)

# Link Eigen librairies
target_link_libraries (RNS_Solver Eigen3::Eigen)

# Add an executable
add_executable (RNS_Solver_restart Load_Cerial.cpp)

# Add an executable

# Link Eigen librairies
target_link_libraries (RNS_Solver_restart Eigen3::Eigen)

# add_executable (RNS_Solver_test Block_Triagonal_Matrix_Inverse_Test.cpp)
# target_link_libraries (RNS_Solver_test Eigen3::Eigen)
