# Require same CMake as MFEM
cmake_minimum_required(VERSION 3.10)

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define our project
project(jots VERSION 1.0
	     LANGUAGES CXX)

# Get preCICE
find_package(PRECICE REQUIRED CONFIG NAMES precice)
message(STATUS "Found preCICE config: TRUE")

# Get MFEM
find_package(MFEM REQUIRED CONFIG)
message(STATUS "Found MFEM config: TRUE")

# Pass the version to the .hpp.in file --> .hpp file in build directory
configure_file(global_vars.hpp.in  global_vars.hpp)


# Add src directory (creates jotsLib)
add_subdirectory(src)

# Create executable
add_executable(jots jots_main.cpp)
target_link_libraries(jots PRIVATE precice::precice)
target_link_libraries(jots PRIVATE mfem) # See the MFEMConfig.cmake file as this is name of library
target_link_libraries(jots PRIVATE jotsLib)

# Include the include file made above from the build directory
# Also include jotsLib files (replicating Step2 tutorial for cmake)
target_include_directories(jots PUBLIC "${PROJECT_BINARY_DIRECTORY}"
					src)
