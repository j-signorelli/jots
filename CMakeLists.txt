# Require same CMake as MFEM
cmake_minimum_required(VERSION 3.10)

# Define our project
project(jots VERSION 1.0)

# Create executable
add_executable(jots src/jots.cpp)

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Pass the version to the .hpp.in file --> .hpp file in build directory
configure_file(include/global_vars.hpp.in  global_vars.hpp)

# Include this include file from the build directory
target_include_directories(jots PUBLIC "${PROJECT_BINARY_DIRECTORY}")

