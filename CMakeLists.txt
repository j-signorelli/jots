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
# if (DEFINED PRECICE_FOUND)
#	add_compile_definitions(PRECICE_ENABLED)
message(STATUS "Found preCICE config: TRUE")
#endif()

# Get MFEM
find_package(MFEM REQUIRED CONFIG)
message(STATUS "Found MFEM config: TRUE")

# Get Boost headers
find_package(Boost REQUIRED CONFIG COMPONENTS headers)
message(STATUS "Found Boost headers: TRUE")

# Add src directory (creates jotsLib)
add_subdirectory(src)

# Create executable
add_executable(jots jots_main.cpp)
target_link_libraries(jots PRIVATE precice::precice)
target_link_libraries(jots PRIVATE mfem) # See the MFEMConfig.cmake file as this is name of library
target_link_libraries(jots PRIVATE jotsLib)

# Also include jotsLib files (replicating Step2 tutorial for cmake)
target_include_directories(jots PUBLIC src)
